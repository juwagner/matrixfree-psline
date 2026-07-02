# ------------------------------------------------------------------------------
# Test suite for high-level additive P-spline matrix-free operators
# ------------------------------------------------------------------------------

rm(list=ls())

library(rTensor)
library(Matrix)

source("src/pspline/pspline_matrices.R")
source("src/pspline_additive/pspline_operations_additive.R")

# ------------------------------------------------------------------------------
# Generate test data
# ------------------------------------------------------------------------------

source("tests/pspline_additive/generate_test_data_additive.R")

mvp_ref <- function(A, x) {
  as.vector(A %*% x)
}

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

# Phi
mvp_Phi_mf <- mvp_Phi_terms(PhiT_terms, alpha_terms)
mvp_Phi_ref <- mvp_ref(Phi, alpha)
cat("MVP with Phi correct:", 
    all.equal(mvp_Phi_mf, mvp_Phi_ref, tol=1e-10), "\n")

# PhiT
mvp_PhiT_mf <- unlist(mvp_PhiT_terms(PhiT_terms, y))
mvp_PhiT_ref <- mvp_ref(PhiT_full, y)
cat("MVP with Phi^T correct:", 
    all.equal(mvp_PhiT_mf, mvp_PhiT_ref, tol=1e-10), "\n")

# mvp_PhiTPhi
mvp_PhiTPhi_mf <- unlist(mvp_PhiTPhi_terms(PhiT_terms, alpha_terms))
mvp_PhiTPhi_ref <- mvp_ref(PhiTPhi_full, alpha)
cat("MVP with Phi^T %*% Phi correct:", 
    all.equal(mvp_PhiTPhi_mf, mvp_PhiTPhi_ref, tol=1e-10), "\n")

# mvp_lambda_Lambda
mvp_Lambda_mf <- unlist(
  mvp_lambda_Lambda_terms(L_terms, lambda_vec, alpha_terms)
)
mvp_Lambda_ref <- mvp_ref(lambda_Lambda_full, alpha)
cat("MVP with (weighted) Lambda correct:", 
    all.equal(mvp_Lambda_mf, mvp_Lambda_ref, tol=1e-10), "\n")

# mvp_A_lambda
mvp_A_lambda_mf <- unlist(
  mvp_A_lambda_terms(PhiT_terms, L_terms, lambda_vec, alpha_terms)
)
mvp_A_lambda_ref <- mvp_ref(A_lambda_full, alpha)
cat("MVP with A_lambda correct:", 
    all.equal(mvp_A_lambda_mf, mvp_A_lambda_ref, tol=1e-10), "\n")
