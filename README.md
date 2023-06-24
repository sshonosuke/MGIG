# Gibbs Sampler for Matrix Generalized Inverse Gaussian Distributions

This repository provides R code implementing Gibbs sampler for matrix generalized inverse Gaussian (MGIG) distributions, as proposed by the following paper.

Hamura, Y., Irie, K. and Sugasawa, S. (2023). Gibbs sampler for matrix generalized inverse Gaussian distribution. *arXiv:2302.09707*.

The repository includes the following files.

- `MGIG-MCMC-function.R` : Implementation of Gibbs sampler and existing sampling algorithms to generate multiple random samples from MGIG
- `MGIG-sampler.R` : Implementation of sampling algorithms for MGIG as one-step update (This can be implemented as a part of MCMC algorithm)
- `PGGM-function.R`: Implementing MCMC algorithm for partial Gaussian Graphical models (PGGM)
- `MVST-function.R`: Implementing MCMC algorithm for matrix-variate skew-$t$ (MVST) distributions 
- `Example-generation.R`: Comparison of sampling algorithms for MGIG
- `Example-PGGM.R`: One-shot simulation study with PGGM
- `Example-MVST.R`: Example of fitting MVST distributions to landsat satellite data
- `sat.trn`: Landsat satellite data
