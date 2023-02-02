# SCTMM

## Getting started

This repository hosts the code for SCTMM, a scalable approach for continuous time Markov models, as described in "A scalable approach for continuous time Markov models with covariates" paper.

## About CTMM

Existing methods for fitting continuous time Markov models (SCTMM) in the presence of covariates suffer from scalability issues due to high computational cost of matrix exponentials calculated for each observation. In this paper we propose an optimization technique for SCTMM which uses a stochastic gradient descent algorithm combined with differentiation of the matrix exponential using a Pad\'e approximation which we show then makes fitting large scale data feasible. We present two novel methods for computing standard errors, one using Pad\'e expansion, and the other using power series expansion of the matrix exponential. Through simulations we find improved performance relative to existing SCTMM methods, and we demonstrate the method on the large-scale multiple sclerosis NO.MS dataset.

There are two directories in this repository:

```
Simulation
SCTMM
```

## Using SCTMM
- The "Simulation/" folder contains the implementation and source code to simulate the dataset needed; two scenarios are generated: first a null case where there is no effect of covariates on the transition rates, and a second case where there is an effect from the covariate impacting transition rates.
- The folder "SCTMM/" contains the source code to reproduce the results of our simulation studies.
   - Specifically, the script `multi_Q_multi_cov.R` run the SCTMM model on simulated data (need to be imported to the same path) and outputs the infinitesimal transition matrix; this involves calculating the gradient of the exponential of a matrix using Pade approximation, and later using a stochastic gradient which is useful for large scale dataset. 
   - The code script `hessian.R` performs the calculation for confidence interval using double usage of Pade approximation. Where needed the calculation in these scripts are parallelized which can be set manually by the user. 

## License
Open source; available for public use.
