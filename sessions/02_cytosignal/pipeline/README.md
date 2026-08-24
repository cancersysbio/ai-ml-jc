# CytoSignal on breast_atera — pipeline

Applies [CytoSignal](https://github.com/welch-lab/cytosignal) (Liu et al.,
*Nature Genetics* 2026) to the `breast_atera` spatial dataset from the
`grafiti` spatial-autoencoder pipeline, clustered by an external breast
cell-type labeling (not grafiti's own motifs — see `CLUSTER_SRC` in the R
script to switch).

This directory holds the supporting scripts behind
[`../implementation.ipynb`](../implementation.ipynb) — the notebook is the
demo/walkthrough; these are what actually produce its input.

## Data

- **AnnData**: a grafiti sweep run's cached inference result (`result.h5ad`)
  — spatial coords in `obsm["spatial"]`, grafiti motif labels in
  `obs["grafiti_spatial"]`, grafiti cell-type cluster ids in
  `obs["grafiti_celltype"]` (e.g. `CT0`..`CT9`)
- **Cell types**: an external pkl mapping those `grafiti_celltype` cluster
  ids to human-readable names (e.g. `{"CT0": "Cancer1_2", ...}`) — **not** a
  per-cell mapping; see `load_celltypes()` in `00_export_breast_atera.py`
- Edit `ADATA_PATH` / `CELLTYPE_PKL` at the top of `00_export_breast_atera.py`
  and `04_export_grafiti_edges.py` (or pass `--adata-path`/`--celltype-pkl`)
  to point at your own copies

## Files, in run order

1. **`00_export_breast_atera.py`** — loads the h5ad + pkl, aligns cell types
   to cells, exports a sparse `matrix.mtx` (genes x cells) + `spatial.csv` +
   `clusters_celltype.csv` (+ `clusters_motif.csv` as a fallback/comparison
   cluster set) for R. Also prints a nearest-neighbor distance diagnostic —
   use it to sanity-check `SCALE_FACTOR` below.

   **Important gotcha**: on this dataset, neither `adata.X` nor
   `adata.layers["counts"]` (despite the name) held raw integer counts — both
   turned out to be `log1p(normalize_total(...))`-processed, all the way back
   through the "raw" source file too. The script applies `expm1()` to undo
   the log step (mathematically exact) before export; this recovers
   library-size-normalized values, not true raw counts, but removes the log
   compression that otherwise makes `removeLowQuality(counts.thresh=...)`
   nonsensical. Check your own data's `adata.layers.keys()` and value ranges
   before assuming this applies to you.
2. **`01_cytosignal_breast_atera.R`** — the CytoSignal pipeline: QC ->
   `inferEpsParams` -> `findNN` -> `imputeLR` -> `inferIntrScore` (the
   expensive permutation-test step) -> `inferSignif` -> `rankIntrSpatialVar`
   -> saves `cytosignal_breast_atera.rds` + top-interaction CSVs + summary
   plots. All parameters are overridable by env var (see the `getenv()`
   block at the top) so the SLURM scripts don't need to edit this file.
   Saves a checkpoint `.rds` right after the expensive scoring step, *before*
   the SPARK-dependent ranking step, and falls back to `"result.hq"`
   significance (instead of `"result.spx"`) if that ranking step fails —
   so a missing/broken optional dependency can't cost you the whole run.
3. **`02_pilot.sbatch`** / **`02_full.sbatch`** — SLURM submission scripts.
   Edit `--partition` for your cluster. **Run the pilot first.** It
   subsamples ~8,000 cells and completes in a bounded time window, which
   tells you the real per-cell cost of `inferIntrScore` on this dataset
   before you commit to a full run — there is no published CytoSignal
   benchmark at hundreds-of-thousands-of-cells scale to size `--time`/`--mem`
   from otherwise. Note that random cell subsampling for a pilot does distort
   the spatial neighbor graph (Delaunay triangulation especially — it
   connects cells that are graph-adjacent in the subsample but weren't truly
   touching in the tissue), so treat pilot *results* as a pipeline sanity
   check only, not a scientific finding. Read the pilot's `.rds` and log,
   adjust `SCALE_FACTOR` / thresholds if `inferSignif()` found 0 significant
   interactions, *then* submit the full run and rescale `--time`/`--mem` off
   the pilot's observed wall-clock and peak memory (`sacct -j <jobid>
   --format=Elapsed,MaxRSS` after it finishes).
4. **`04_export_grafiti_edges.py`** *(optional)* — exports grafiti's own
   learned per-edge interaction scores (`get_edge_scores()`'s
   `edge_interaction_local`/`_global`, and/or attention weights if you want
   them) for comparing CytoSignal's significant LR edges against what
   grafiti's model independently attended to. Only useful if you already
   have a trained grafiti model and its cached `result.h5ad` includes these
   `obsp` matrices (check `adata.obsp.keys()`); otherwise they need
   computing via `model.get_edge_scores(adata)` on a GPU node first.

The results notebook (`../implementation.ipynb`) is the live demo — it loads
a `.rds` produced by step 2/3 above and does not re-run any expensive step.

## One-time R environment setup

CytoSignal is an R package — it needs an R installation, separate from the
Python venv used for the export step. Installing it into its own conda env
(rather than a module-loaded R with a hard-to-control library path) keeps it
alongside your other conda/venv environments:

```bash
conda create -n cytosignal-r -c conda-forge r-base=4.3 r-devtools r-irkernel \
    r-ggplot2 r-knitr r-dplyr r-matrix -y
conda activate cytosignal-r
```
```r
devtools::install_github("welch-lab/cytosignal")
IRkernel::installspec(name = "ir", displayname = "R (cytosignal-r)", user = TRUE)
```

`IRkernel::installspec(user = TRUE)` registers the kernel to
`~/.local/share/jupyter/kernels/`, so any Jupyter server you launch (from any
env) can see and select it — Jupyter itself doesn't need to run inside
`cytosignal-r`.

The `.sbatch` scripts activate this same env by name
(`${CYTOSIGNAL_ENV:-cytosignal-r}`) before running `01_cytosignal_breast_atera.R`
— if you name the env something else, export `CYTOSIGNAL_ENV=<name>` before
submitting, or edit the two `.sbatch` files.

### Gotchas: optional CytoSignal dependencies

`devtools::install_github()` only installs `Imports`, not `Suggests` — but
several `plotSignif()`/`plotCircosNIntr()` code paths this notebook actually
exercises need packages that are only in `Suggests`. Install these too, or
you'll hit them one at a time mid-demo:

```r
install.packages(c("cowplot", "patchwork", "scattermore", "png", "plot3D", "circlize"))
```
- `cowplot` — `plotSignif()`'s multi-panel layout (hard requirement for the
  main results plot, not actually optional in practice)
- `plot3D` — only needed for `plotSignif(..., edge = TRUE)`'s 3D edge plot;
  set `edge = FALSE` to skip it entirely if you don't want the dependency
- `circlize` — only for `plotCircosNIntr()`; also has an unrelated real bug
  in the installed package as of this writing (`countEdges()` uses `%||%`
  unqualified without importing it from `rlang`, which only exists as a base
  R operator from R>=4.4 — a problem if you're on R 4.3.x). If you hit
  `could not find function "%||%"`, don't fight the installed package;
  reimplement the small aggregation function yourself (it's ~20 lines, see
  the notebook for a working version) rather than patching a locked package
  namespace live.
- `png`/`scattermore`/`patchwork` — various plotting/rasterization paths;
  `png` in particular can fail to compile with `cannot find -lz` in a conda
  env missing zlib — fix with `conda install -c conda-forge zlib libpng -y`
  before retrying `install.packages("png")`.
- `rankIntrSpatialVar()` needs `SPARK` (`devtools::install_github('xzhoulab/SPARK')`),
  which as of this writing fails to compile against modern `RcppArmadillo`
  (needs C++14, SPARK's own `Makevars` pins C++11) — fix by patching SPARK's
  `src/Makevars` to `CXX_STD = CXX14` before installing from a local clone,
  since a personal `~/.R/Makevars`/`R_MAKEVARS_USER` override does **not**
  win over a package's own `Makevars`. This is optional — the R script
  already treats a `rankIntrSpatialVar()` failure as non-fatal.

Get Jupyter running with the `ir` kernel visible however your cluster
normally exposes Jupyter (a shared JupyterHub if you have one, or a manual
`jupyter notebook --no-browser` + SSH tunnel from your laptop otherwise);
open `../implementation.ipynb` and pick the "R (cytosignal-r)" kernel.

## Running

```bash
sbatch 02_pilot.sbatch
squeue -u $USER                   # watch it queue/run
tail -f logs/pilot_*.out
```

Once the pilot's `.rds` lands, open the notebook and run through it —
that's your dry run. Submit `02_full.sbatch` once you've adjusted its
resource request from the pilot's numbers (`sacct -j <jobid>
--format=Elapsed,MaxRSS`); swap `RDS_PATH` in the notebook's first code cell
to point at the full output when it's ready.

## Known gaps (not implemented here)

- VeloCytoSignal (needs RNA velocity, not computed for this dataset)
- Multi-sample/patient differential signaling (need per-patient sample IDs
  if your dataset has multiple patients)
- GO enrichment on the signaling-associated DEGs — the notebook stops at the
  significant gene list; wire up `gprofiler2::gost()` or `revigo()` per the
  package vignette if useful
