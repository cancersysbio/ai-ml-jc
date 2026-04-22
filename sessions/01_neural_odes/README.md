# Session 01: Neural Ordinary Differential Equations

- **Date:** 2026-04-20
- **Presenter:** @tulerpetontidae (Artem)
- **Paper:** [Neural Ordinary Differential Equations](https://arxiv.org/abs/1806.07366) (Chen et al., NeurIPS 2018)

## Summary

Instead of specifying a discrete sequence of hidden layers, Neural ODEs parameterize the derivative of the hidden state using a neural network. The output is computed using a black-box ODE solver. This yields continuous-depth models with constant memory cost, adaptive computation, and invertible normalizing flows.

## Key Concepts

- Residual networks as Euler discretization of a continuous transformation
- Adjoint sensitivity method for memory-efficient backpropagation through ODE solvers
- Continuous normalizing flows (CNFs) — generative models without partitioning or ordering constraints
- Latent ODE models for irregularly-sampled time series

## Resources

- [arXiv](https://arxiv.org/abs/1806.07366)
- [Original code (torchdiffeq)](https://github.com/rtqichen/torchdiffeq)

## Implementation

See [`implementation.ipynb`](implementation.ipynb) for a walkthrough.
