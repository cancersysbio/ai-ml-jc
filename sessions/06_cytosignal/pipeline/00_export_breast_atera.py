"""Export breast_atera grafiti result.h5ad + external cell-type pkl to R-readable
files for CytoSignal.

Inputs (edit ADATA_PATH / CELLTYPE_PKL below, or pass --adata-path / --celltype-pkl):
  ADATA_PATH   = grafiti sweep result.h5ad (breast_atera, m12 motif model)
  CELLTYPE_PKL = external cell-type labels for the same run

Outputs (written to OUTPUT_DIR), R/CytoSignal-ready:
  matrix.mtx, barcodes.tsv, features.tsv   -- sparse counts, genes x cells (R convention)
  spatial.csv                              -- columns x, y, one row per cell, same order as barcodes.tsv
  clusters_celltype.csv                    -- column cell_type, from CELLTYPE_PKL (primary clusters input)
  clusters_motif.csv                       -- column grafiti_spatial, from adata.obs (kept for comparison)

Run in a Python env with anndata/scipy (no R needed here):
    python 00_export_breast_atera.py --output-dir ./exported [--subset 8000]
"""
import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from scipy.spatial import cKDTree

ADATA_PATH = "/path/to/your/grafiti_models/breast_atera/<run_name>/result.h5ad"
CELLTYPE_PKL = "/path/to/your/breast_cts_labels/<run_name>.pkl"


def load_celltypes(pkl_path: str, adata) -> pd.Series:
    """Map adata.obs['grafiti_celltype'] (cluster ids, e.g. 'CT0'..'CT9')
    through the external pkl (cluster id -> human-readable name) to get a
    per-cell cell-type label.

    Expected pkl shape: a small dict like {'CT0': 'Cancer1_2', 'CT1':
    'Cancer1_1', ...}, one entry per adata.obs['grafiti_celltype'] cluster id
    -- NOT a per-cell mapping.
    """
    with open(pkl_path, "rb") as fh:
        name_by_cluster = pickle.load(fh)
    if not isinstance(name_by_cluster, dict):
        raise TypeError(
            f"Expected a dict mapping grafiti_celltype id -> name, got {type(name_by_cluster)}"
        )

    if "grafiti_celltype" not in adata.obs:
        raise KeyError("adata.obs['grafiti_celltype'] not found -- cannot map the "
                        "cluster-id -> name pkl onto cells")
    cluster_ids = adata.obs["grafiti_celltype"].astype(str)

    missing_keys = sorted(set(cluster_ids.unique()) - set(name_by_cluster))
    if missing_keys:
        raise ValueError(
            f"grafiti_celltype has cluster ids not present in the pkl: {missing_keys} "
            f"(pkl has: {sorted(name_by_cluster)})"
        )

    s = cluster_ids.map(name_by_cluster)
    s.index = adata.obs_names
    print(f"[celltype pkl] mapped {len(name_by_cluster)} clusters -> "
          f"{s.nunique()} cell-type names across {len(s)} cells")
    print(f"[celltype pkl] counts:\n{s.value_counts().to_string()}")
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata-path", default=ADATA_PATH)
    ap.add_argument("--celltype-pkl", default=CELLTYPE_PKL)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--subset", type=int, default=None,
                     help="Randomly sample N cells for a pilot run (default: full dataset)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    import anndata as ad

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.adata_path} ...")
    adata = ad.read_h5ad(args.adata_path)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    if "spatial" not in adata.obsm:
        sys.exit("adata.obsm['spatial'] not found -- this script assumes grafiti's "
                  "squidpy-convention coords; check adata.obsm.keys()")
    coords = np.asarray(adata.obsm["spatial"])[:, :2]

    celltype = load_celltypes(args.celltype_pkl, adata)

    if "grafiti_spatial" in adata.obs:
        motif = adata.obs["grafiti_spatial"].astype(str)
    else:
        print("WARNING: adata.obs['grafiti_spatial'] not found -- skipping motif export")
        motif = None

    keep = celltype.notna().to_numpy()
    if not keep.all():
        print(f"Dropping {(~keep).sum()} cells with missing cell-type label")
        adata = adata[keep].copy()
        coords = coords[keep]
        celltype = celltype[keep]
        if motif is not None:
            motif = motif[keep]

    if args.subset is not None and args.subset < adata.n_obs:
        rng = np.random.default_rng(args.seed)
        idx = np.sort(rng.choice(adata.n_obs, size=args.subset, replace=False))
        adata = adata[idx].copy()
        coords = coords[idx]
        celltype = celltype.iloc[idx]
        if motif is not None:
            motif = motif.iloc[idx]
        print(f"Subsampled to {adata.n_obs} cells (--subset {args.subset}, seed={args.seed})")

    # --- counts: genes x cells (R convention), sparse ---
    # grafiti convention: raw counts live in adata.layers["counts"] (used for
    # HVG/scVI upstream); adata.X may hold normalized/log expression instead.
    # CytoSignal's removeLowQuality(counts.thresh=...) needs raw counts.
    # Confirmed on this dataset: adata.X, adata.layers["counts"], and even the
    # "raw" source h5ad all hold the SAME log1p(normalize_total(...)) values,
    # subset to these 4000 genes -- there is no true raw-count matrix cached
    # anywhere in this chain (obs["transcript_counts"]/["total_counts"] are the
    # only place the real per-cell raw totals survive, as scalars, not per-gene).
    # expm1() exactly undoes log1p (expm1(log1p(x)) == x), recovering the
    # normalize_total'd (but unlogged) values -- NOT raw integer counts, but at
    # least removes the log compression that made counts.thresh nonsensical.
    if "counts" in adata.layers:
        X = adata.layers["counts"]
    else:
        print("WARNING: adata.layers['counts'] not found -- falling back to adata.X")
        X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocoo()
    X.data = np.expm1(X.data)
    print("Applied expm1() to undo log1p (data is library-size-normalized, not raw counts)")
    # AnnData is cells x genes; CytoSignal wants genes x cells
    counts_gxc = sp.coo_matrix((X.data, (X.col, X.row)), shape=(adata.n_vars, adata.n_obs)).tocsr()

    print(f"Writing matrix.mtx ({counts_gxc.shape[0]} genes x {counts_gxc.shape[1]} cells) ...")
    sio.mmwrite(str(out / "matrix.mtx"), counts_gxc)
    pd.Series(adata.obs_names).to_csv(out / "barcodes.tsv", index=False, header=False)
    pd.Series(adata.var_names).to_csv(out / "features.tsv", index=False, header=False)

    pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=adata.obs_names).to_csv(
        out / "spatial.csv"
    )
    pd.DataFrame({"cell_type": celltype.to_numpy()}, index=adata.obs_names).to_csv(
        out / "clusters_celltype.csv"
    )
    if motif is not None:
        pd.DataFrame({"grafiti_spatial": motif.to_numpy()}, index=adata.obs_names).to_csv(
            out / "clusters_motif.csv"
        )

    # --- diagnostics: nearest-neighbor distance, to sanity-check scale.factor in R ---
    tree = cKDTree(coords)
    nn_dist, _ = tree.query(coords, k=2)  # k=1 is self (dist 0)
    nn_dist = nn_dist[:, 1]
    print("\n=== Summary ===")
    print(f"n cells: {adata.n_obs}, n genes: {adata.n_vars}")
    print(f"coordinate range: x [{coords[:,0].min():.1f}, {coords[:,0].max():.1f}], "
          f"y [{coords[:,1].min():.1f}, {coords[:,1].max():.1f}]")
    print(f"nearest-neighbor distance (raw coord units): "
          f"median={np.median(nn_dist):.2f}, p5={np.percentile(nn_dist,5):.2f}, "
          f"p95={np.percentile(nn_dist,95):.2f}")
    print("  -> typical cell-cell spacing in tissue is ~10-20 um; compare against the "
          "median above to sanity-check SCALE_FACTOR in 01_cytosignal_breast_atera.R "
          "(SCALE_FACTOR = um_per_unit; if coords are already in um, SCALE_FACTOR = 1.0)")
    print(f"unique cell types: {sorted(celltype.astype(str).unique())}")
    if motif is not None:
        print(f"unique grafiti motifs: {sorted(motif.astype(str).unique())}")
    print(f"\nWrote exported files to {out}/")


if __name__ == "__main__":
    main()
