# Session 02: Nicheformer - Transformers for Single Cell Transcriptomics

- **Date:** 2026-05-13
- **Presenter:** @bsimon11 (Brennan)
- **Paper:** [Nicheformer: a foundation model for single-cell and spatial omics]([https://arxiv.org/abs/1806.07366](http://nature.com/articles/s41592-025-02814-z)) (Tejada-Lapuerta et al., Nat Methods 2025)

## Summary

Nicheformer is a pre-trained transformer-based foundation model designed to embed single-cell RNA sequencing data from diverse tissue samples, disease contexts, and technologies within a common latent space. 

## Key Concepts

- Nicheformer embeds single-cell RNAseq data into a common latent space
- It is pretrained on many technologies as well as spatial and non-spatial, but does not include CosMx 6k or WT (it does include CosMx 1K)
- In our data, Nicheformer embeddings do not cluster by cell type, possibly driven by sample quality
- This type of model is promising for single-cell analysis, but perhaps doesn't outperform sample-specific UMAP yet (for CosMx 6k data)

## Resources

- [paper]([https://arxiv.org/abs/1806.07366](https://www.nature.com/articles/s41592-025-02814-z))
- [Original code (nicheformer)]([https://github.com/rtqichen/torchdiffeq](https://github.com/theislab/nicheformer))

## Implementation

See get_embeddingsTMA1.py and ViewEmbsClean.ipynb for a walkthrough.
