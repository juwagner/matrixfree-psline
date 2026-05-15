# Matrix-free Penalized Spline Smoothing with Multiple Covariates

This repository contains the code to the paper:

**Matrix-free Penalized Spline Smoothing with Multiple Covariates**  
Julian Wagner, Göran Kauermann, Ralf Münnich  
arXiv: https://arxiv.org/abs/2101.06034

Note that it is researches based code in order to comprehend the methods provided within the paper.
It is no production ready code and should be used with caution.

## Overview

This repository implements the matrix-free penalized spline smoothing approach described in the paper.  
The main contributions reproduced here are:

- Efficient spline smoothing over multiple covariates without constructing large tensor-product bases
- Matrix-free computation of penalized spline estimators
- Mixed-model formulation enabling computationally efficient smoothing parameter estimation
- Extensions towards generalized regression models
- Example applications, including high-dimensional spatial and satellite-based data


## Open work

- Make C++ methods independent of P
- Improve performance of C++ methods
- Add correct diagonal precontionier for generalized models
- Make PCG methods numerically stable
- Add stopping criterion for fixpoint iteration