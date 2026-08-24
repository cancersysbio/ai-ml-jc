"""Export grafiti's edge_interaction_local (from get_edge_scores()) for
comparison against CytoSignal's sender/receiver significance calls.

Uses the local (softmax-within-target-neighborhood) score, per grafiti's own
docstring the one meant for "attribute influence within a cell's
neighborhood" -- the most directly comparable quantity to CytoSignal's
per-edge signaling contribution. edge_interaction_global (uncalibrated MLP
output, comparable across the whole tissue) is exported too, in case you want
a tissue-wide ranking instead of a per-neighborhood one.

Run in a Python env with anndata/scipy (the same one used for
00_export_breast_atera.py, or any env with grafiti's dependencies):
    python 04_export_grafiti_edges.py --output-dir ./exported_full
"""
import argparse
from pathlib import Path

import anndata as ad
import scipy.io as sio

ADATA_PATH = "/path/to/your/grafiti_models/breast_atera/<run_name>/result.h5ad"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata-path", default=ADATA_PATH)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.adata_path} ...")
    adata = ad.read_h5ad(args.adata_path)

    for key in ["edge_interaction_local", "edge_interaction_global"]:
        if key not in adata.obsp:
            print(f"WARNING: {key} not found in adata.obsp, skipping")
            continue
        mtx = adata.obsp[key]
        print(f"{key}: {mtx.shape}, {mtx.nnz} nonzero edges")
        sio.mmwrite(str(out / f"grafiti_{key}.mtx"), mtx)

    # barcodes for this matrix's row/col order -- should match exported_full/barcodes.tsv
    # already written by 00_export_breast_atera.py off the same adata, but write our own
    # copy here too so this script is usable standalone.
    with open(out / "grafiti_edge_barcodes.tsv", "w") as f:
        f.write("\n".join(adata.obs_names))

    print(f"Wrote grafiti edge matrices to {out}/")


if __name__ == "__main__":
    main()
