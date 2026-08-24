# Session 06: CytoSignal

- **Date:** 2026-08-24
- **Presenter:** Maryam Pourmaleki
- **Paper:** [CytoSignal Detects Locations and Dynamics of Ligand-Receptor Signaling at Cellular Resolution from Spatial Transcriptomic Data](https://www.nature.com/articles/s41588-026-02624-9) (Liu et al., *Nature Genetics* 2026)

## Summary

CytoSignal infers spatially-resolved, cell-cell ligand-receptor signaling
from spatial transcriptomics at single-cell resolution. It separately models
two signaling modes with different spatial neighborhoods: diffusion-dependent
interactions (Gaussian epsilon-ball, for cytokines/chemokines/growth factors
reaching cells at a distance) and contact-dependent interactions (Delaunay
triangulation, for juxtacrine signaling between physically touching cells).
For each candidate ligand-receptor pair (from a CellPhoneDB-derived
database), it computes a spatially-resolved co-expression score and tests
significance via a spatial permutation test with spatial FDR correction.

## Key Concepts

- Epsilon-ball vs. Delaunay-triangulation neighborhoods for diffusion- vs.
  contact-dependent signaling
- Spatial permutation testing + spatial FDR correction for calling
  significant signaling locations
- SPARK-X spatial-variability ranking of significant interactions
- Signaling-associated differentially expressed genes via elastic-net
  regression on the interaction test statistic

## Resources

- [Paper (Nature Genetics)](https://www.nature.com/articles/s41588-026-02624-9)
- [CytoSignal GitHub](https://github.com/welch-lab/cytosignal)
- [Slides](https://docs.google.com/presentation/d/1O0ouGGh-Aty_DyArXxZWSMcRJXdG3iKNs5D55qXKOCA/edit?slide=id.p#slide=id.p)

## Implementation

Applied to a real breast cancer spatial dataset (`breast_atera`, ~170K
cells; originally [10x Genomics' Atera WTA FFPE Human Breast Cancer
dataset](https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer)),
clustered using an external cell-type labeling of a
[grafiti](https://github.com/nceglia/grafiti) spatial-autoencoder model run
on the same data. See [`implementation.ipynb`](implementation.ipynb) for the
live results walkthrough, and [`pipeline/`](pipeline/) for the scripts that
produce its input (export from AnnData, the CytoSignal R pipeline, and SLURM
submission scripts) — see [`pipeline/README.md`](pipeline/README.md) for the
full pipeline write-up, including gotchas hit along the way (raw-counts
availability, optional-dependency installation issues, and a couple of real
bugs found in CytoSignal's own `Suggests`-gated code paths).
