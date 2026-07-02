# ------------------------------------------------------------------------------
# Test suite for PCG solver + diagonal matrix-free operators for additive P-splines
# ------------------------------------------------------------------------------

rm(list=ls())

library(rTensor)
library(Matrix)

source("src/pspline/pspline_matrices.R")
source("src/pspline_additive/pcg_solver_additive.R")

# ------------------------------------------------------------------------------
# Generate test data
# ------------------------------------------------------------------------------

source("tests/pspline_additive/generate_test_data_additive.R")

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

b_terms = mvp_PhiT_terms(PhiT_terms = PhiT_terms, x = y)

alpha_pcg <- solve_pcg_terms(
  PhiT_terms = PhiT_terms, 
  L_terms = L_terms, 
  lambda_vec = lambda_vec, 
  b_terms = b_terms,
  tol = 1e-12,
  verbose=TRUE
)

b <- PhiT_full %*% y
alpha_full <- as.vector(solve(A_lambda_full, b))

y_hat_pcg <- Phi %*% unlist(alpha_pcg)
y_hat_full <- Phi %*% alpha_full

cat("PCG solver predictions correct: ", 
    all.equal(y_hat_pcg, y_hat_full, tol=1e-10),  "\n")
