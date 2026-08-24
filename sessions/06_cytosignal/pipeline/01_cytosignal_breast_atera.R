# CytoSignal on breast_atera (grafiti sweep run
# breast_atera_nosplit_nobatch__large_lr1e-4_m12_eps0.01_w1000), clustered by
# an external breast_cts cell-type labels pkl for the same run.
#
# Input: the directory written by 00_export_breast_atera.py (matrix.mtx,
# barcodes.tsv, features.tsv, spatial.csv, clusters_celltype.csv[,
# clusters_motif.csv]).
#
# Usage:
#   Rscript 01_cytosignal_breast_atera.R  (reads the === USER-CONFIGURABLE
#   PARAMETERS === block below; edit DATA_DIR / PILOT before submitting)

# === USER-CONFIGURABLE PARAMETERS ===
# Each can be overridden by an env var of the same name (set by the SLURM
# scripts) so 02_pilot.sbatch / 02_full.sbatch don't need to edit this
# file -- e.g. `DATA_DIR=./exported_pilot Rscript 01_cytosignal_breast_atera.R`
getenv <- function(name, default) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v)) return(default)
  if (is.numeric(default)) return(as.numeric(v))
  if (is.logical(default)) return(as.logical(v))
  v
}

DATA_DIR     <- getenv("DATA_DIR", "./exported")          # output of 00_export_breast_atera.py
OUTPUT_DIR   <- getenv("OUTPUT_DIR", "./cytosignal_output")
CLUSTER_SRC  <- getenv("CLUSTER_SRC", "celltype")          # "celltype" (default) or "motif"
SCALE_FACTOR <- getenv("SCALE_FACTOR", 1.0)                # um per coordinate unit -- VERIFY
                                                            # against the nearest-neighbor
                                                            # diagnostic from 00_export_breast_atera.py
PERM_SIZE    <- getenv("PERM_SIZE", 1e5)                   # inferIntrScore's perm.size (package default)
P_THRESH     <- getenv("P_THRESH", 0.05)
READS_THRESH <- getenv("READS_THRESH", 100)
SIG_THRESH   <- getenv("SIG_THRESH", 100)
SEED         <- getenv("SEED", 42)
N_CORES      <- getenv("N_CORES", 8)                       # match --cpus-per-task
COUNTS_THRESH <- getenv("COUNTS_THRESH", 10)                # NOTE: the exported matrix is
                                                              # expm1(log1p(normalize_total(...)))
                                                              # -- library-size-normalized, NOT raw
                                                              # counts -- so this is a near-no-op
                                                              # QC floor, not a real raw-count cutoff.
                                                              # See 00_export_breast_atera.py.

# === Setup ===
suppressPackageStartupMessages({
  library(cytosignal)
  library(Matrix)
})

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
cat(sprintf("[%s] Starting breast_atera CytoSignal run (DATA_DIR=%s)\n", t0, DATA_DIR))

# === 1. Load exported data ===
cat("Loading counts matrix...\n")
dge <- as(Matrix::readMM(file.path(DATA_DIR, "matrix.mtx")), "CsparseMatrix")
barcodes <- readLines(file.path(DATA_DIR, "barcodes.tsv"))
features <- readLines(file.path(DATA_DIR, "features.tsv"))
rownames(dge) <- features
colnames(dge) <- barcodes

spatial <- read.csv(file.path(DATA_DIR, "spatial.csv"), row.names = 1)
spatial <- as.matrix(spatial[barcodes, c("x", "y")])
colnames(spatial) <- c("x", "y")

cluster_file <- if (CLUSTER_SRC == "celltype") "clusters_celltype.csv" else "clusters_motif.csv"
cluster_col  <- if (CLUSTER_SRC == "celltype") "cell_type" else "grafiti_spatial"
clusters_df <- read.csv(file.path(DATA_DIR, cluster_file), row.names = 1)
cluster <- factor(clusters_df[barcodes, cluster_col])
names(cluster) <- barcodes

cat(sprintf("Loaded %d genes x %d cells, %d cluster levels (%s): %s\n",
            nrow(dge), ncol(dge), nlevels(cluster), CLUSTER_SRC,
            paste(levels(cluster), collapse = ", ")))

# === 2. Create CytoSignal object + LR database ===
cs <- createCytoSignal(raw.data = dge, cells.loc = spatial, clusters = cluster)
cs <- addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)

# === 3. QC ===
cs <- removeLowQuality(cs, counts.thresh = COUNTS_THRESH)
cs <- changeUniprot(cs)

# === 4. Spatial neighborhoods ===
cs <- inferEpsParams(cs, scale.factor = SCALE_FACTOR)
print(cs)  # inspect inferred sigma before committing to the full run
cs <- findNN(cs)

# === 5. Impute L/R expression in neighborhoods ===
cs <- imputeLR(cs)

# === 6. Score + permutation test (expensive step) ===
t.score.start <- Sys.time()
cat(sprintf("[%s] Starting inferIntrScore (perm.size=%.0e, numCores=%d)...\n",
            t.score.start, PERM_SIZE, N_CORES))
set.seed(SEED)
cs <- inferIntrScore(cs, perm.size = PERM_SIZE, numCores = N_CORES)
t.score.end <- Sys.time()
cat(sprintf("[%s] inferIntrScore done, elapsed: %s\n", t.score.end,
            format(t.score.end - t.score.start)))

# === 7. Significance calling ===
cs <- inferSignif(cs, p.thresh = P_THRESH, reads.thresh = READS_THRESH, sig.thresh = SIG_THRESH)

# --- Checkpoint: save immediately after the expensive steps (inferIntrScore +
# inferSignif), BEFORE anything that can still fail (e.g. rankIntrSpatialVar
# needing the SPARK package). Losing 30-60+ min of compute to a missing
# package on the last step is a real failure mode we hit once already. ---
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(cs, file.path(OUTPUT_DIR, "cytosignal_breast_atera_checkpoint.rds"))
cat(sprintf("Checkpoint saved: %s\n", file.path(OUTPUT_DIR, "cytosignal_breast_atera_checkpoint.rds")))

# === 8. Spatial variability ranking (creates the "result.spx" significance level) ===
# Non-fatal: if SPARK isn't installed (or this step fails for any reason), fall
# back to "result.hq" for everything downstream instead of losing the run.
SIGNIF_USE <- "result.spx"
spx_ok <- tryCatch({
  cs <- rankIntrSpatialVar(cs)
  TRUE
}, error = function(e) {
  cat(sprintf("\n*** rankIntrSpatialVar() failed: %s ***\n", conditionMessage(e)))
  cat("Continuing with signif.use = 'result.hq' instead of 'result.spx' ",
      "(no spatial-variability ranking, but interactions/plots are unaffected ",
      "other than sort order).\n\n", sep = "")
  FALSE
})
if (!spx_ok) SIGNIF_USE <- "result.hq"

n_diff <- length(showIntr(cs, slot.use = "diffusion-Raw_smooth", signif.use = SIGNIF_USE))
n_cont <- length(showIntr(cs, slot.use = "contact-Raw_smooth", signif.use = SIGNIF_USE))
if (n_diff == 0 && n_cont == 0) {
  cat("\n*** DIAGNOSTIC: 0 significant interactions after inferSignif(). ***\n")
  cat("Try lowering READS_THRESH / SIG_THRESH, or check SCALE_FACTOR against the\n")
  cat("nearest-neighbor distance diagnostic from 00_export_breast_atera.py --\n")
  cat("a wrong scale factor makes the epsilon-ball radius nonsensical.\n\n")
}

# === 9. Save full object ===
saveRDS(cs, file.path(OUTPUT_DIR, "cytosignal_breast_atera.rds"))
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "cytosignal_breast_atera.rds")))

# === 10. Export top results tables ===
top_diff <- showIntr(cs, slot.use = "diffusion-Raw_smooth", signif.use = SIGNIF_USE, return.name = TRUE)
top_cont <- showIntr(cs, slot.use = "contact-Raw_smooth", signif.use = SIGNIF_USE, return.name = TRUE)
write.csv(data.frame(id = names(top_diff), name = top_diff)[seq_len(min(20, length(top_diff))), ],
          file.path(OUTPUT_DIR, "top20_diffusion.csv"), row.names = FALSE)
write.csv(data.frame(id = names(top_cont), name = top_cont)[seq_len(min(20, length(top_cont))), ],
          file.path(OUTPUT_DIR, "top20_contact.csv"), row.names = FALSE)

# === 11. Summary plots ===
if (length(top_diff) > 0) {
  plotSignif(cs, intr = seq_len(min(3, length(top_diff))), slot.use = "diffusion-Raw_smooth",
             signif.use = SIGNIF_USE, plot_dir = file.path(OUTPUT_DIR, "plots_diffusion"))
}
if (length(top_cont) > 0) {
  plotSignif(cs, intr = seq_len(min(3, length(top_cont))), slot.use = "contact-Raw_smooth",
             signif.use = SIGNIF_USE, plot_dir = file.path(OUTPUT_DIR, "plots_contact"))
}
# Non-fatal: plotCircosNIntr() needs the 'circlize' package, which is only in
# cytosignal's Suggests (not installed by a default devtools::install_github),
# so don't let a missing optional package take down a run that already saved
# successfully. Also note the real arg names are lrscore.use=/intr.use=, not
# slot.use=/intr= (those belong to plotSignif, not this function).
if (length(top_diff) > 0) {
  png(file.path(OUTPUT_DIR, "circos_top_diffusion.png"), width = 1600, height = 1600, res = 200)
  tryCatch({
    plotCircosNIntr(cs, intr.use = names(top_diff)[1], lrscore.use = "diffusion-Raw_smooth")
  }, error = function(e) {
    cat(sprintf("\nCircos plot skipped (%s). Install circlize with ",
                conditionMessage(e)))
    cat("install.packages('circlize') if you want it for the demo.\n\n")
  }, finally = dev.off())
}

t1 <- Sys.time()
cat(sprintf("\n[%s] Done. Total elapsed: %s\n", t1, format(t1 - t0)))
sessionInfo()
