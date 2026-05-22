# ------------------------------------------------------------------------------
# P-spline models for the LUCAS dataset
# ------------------------------------------------------------------------------

rm(list=ls())

library(Rcpp)
library(tictoc)

sourceCpp("src/base/matrix_free_operations.cpp")
source("src/pspline/pspline_matrices.R")
source("src/pspline/parameter_estimation.R")
source("src/generalized_pspline/generalized_pspline_operations.R")
source("src/generalized_pspline/generalized_parameter_estimation.R")

# ------------------------------------------------------------------------------
# load data
source("src/utils/load_lucas_data.R")

X_input <- X           # U
main <- "y = s(x)"     # s(u)

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
# Solve for α using a fixed λ (fixpoint iteration)

lambda <- 0.1
n_iter <- 3

tic("Fixpoint iteration for generalized p-spline")

alpha <- fixpoint_w_alpha(
  n_iter=n_iter,
  PhiT_list=PhiT_list, 
  L_list=L_list,
  lambda=lambda, 
  pcg_tol=10^(-2),
  pcg_verbose=FALSE
)

toc()

# ------------------------------------------------------------------------------
# Solve for α and λ

M <- 3
V_rad <- rademacher_matrix(K, M, seed=42)
n_iter <- 2

lambda <- 0.1
alpha <- NULL

tic("Fixpoint iteration for generalized p-spline (α and λ)")

for (i in 1:2) {
  cat("Iteration: ", i, "\n")
  alpha <- fixpoint_w_alpha(
    n_iter=n_iter,
    PhiT_list=PhiT_list, 
    L_list=L_list,
    lambda=lambda, 
    #alpha_init=alpha,
    pcg_tol=10^(-2)
  )
  cat("Solved for alpha \n")
  lambda <- estimate_w_lambda(
    n_iter=n_iter,
    PhiT_list=PhiT_list, 
    L_list=L_list,
    alpha=alpha,
    lambda=lambda,
    V_rad=V_rad,
    pcg_tol=10^(-2)
  )
  cat("Current lambda: ", lambda, "\n")
}

toc()

