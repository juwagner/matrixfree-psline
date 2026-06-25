# ------------------------------------------------------------------------------
# P-spline models for the LUCAS dataset
# ------------------------------------------------------------------------------

rm(list=ls())

library(Rcpp)
library(tictoc)

sourceCpp("src/base/matrix_free_operations.cpp")
source("src/pspline/pspline_matrices.R")
source("src/pspline_generalized/pspline_operations_generalized.R")
source("src/pspline_generalized/parameter_estimation_generalized.R")

# ------------------------------------------------------------------------------
# load data
source("src/utils/load_lucas_data.R")

X_input <- X           # U

# ------------------------------------------------------------------------------
# P-spline setup

P     <- ncol(X_input)        # dimension of tensor product
m     <- rep(36, P)           # interior knots per dimension
q     <- rep(3,  P)           # spline degree per dimension
l     <- rep(2,  P)           # difference penalty order
J_vec <- m + q + 1            # basis sizes
K     <- prod(J_vec)          # total number of coefficients


PhiT_list <- lapply(
  1:P, function(p) build_univarate_bspline_basis_T(X_input[,p], m[p], q[p])
)

L_list <- lapply(1:P, function(p) build_penalty_difference(J=J_vec[p], l=l[p]))

# ------------------------------------------------------------------------------
# Estimate α using a fixed λ

lambda <- 0.1
n_iter <- 3

tic("Iteration to estimate alpha for generalized p-spline")

alpha <- estimate_alpha_generalized(
  n_iter=n_iter,
  PhiT_list=PhiT_list, 
  L_list=L_list,
  lambda=lambda, 
  pcg_tol=10^(-4),
  pcg_verbose=FALSE
)

toc()

# ------------------------------------------------------------------------------
# Solve for α and λ

V_rad <- rademacher_matrix(K, M=3, seed=42)
lambda_init <- 0.1

tic("Fixpoint iteration for generalized p-spline (α and λ)")

result <- fit_pspline_generalized(
    n_iter=2,
    n_iter_alpha=2,
    n_iter_lambda=2,
    PhiT_list=PhiT_list,
    L_list=L_list,
    lambda=lambda_init,
    V_rad,
    pcg_tol=10^(-2)
)

alpha <- result$alpha
lambda <- result$lambda

toc()

