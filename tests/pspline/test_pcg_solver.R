# ------------------------------------------------------------------------------
# Test suite for PCG solver + diagonal matrix-free operators for P-splines
# ------------------------------------------------------------------------------

rm(list=ls())

library(microbenchmark)

source("src/pspline/pspline_matrices.R")
source("src/pspline/pcg_solver.R")

# ------------------------------------------------------------------------------
# Generate test data
# ------------------------------------------------------------------------------

source("tests/pspline/generate_test_data.R")

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

# diag_Lambda
diag_Lambda_mf <- get_diag_Lambda(L_list=L_list)
diag_Lambda <- diag(Lambda)
cat("Diagonal of Lambda correct:", 
    all.equal(diag_Lambda_mf, diag_Lambda, tol=1e-10), "\n")

# diag_PhiTPhi
diag_PhiTPhi_mf <- get_diag_PhiTPhi(PhiT_list=PhiT_list)
diag_PhiTPhi <- diag(PhiTPhi)
cat("Diagonal of PhiTPhi correct:", 
    all.equal(diag_PhiTPhi_mf, diag_PhiTPhi, tol=1e-10), "\n")

# PCG solver
b <- as.vector(PhiT %*% y)
alpha_pcg  <- solve_pcg(PhiT_list, L_list, lambda, b, tol = 1e-12, verbose=TRUE)
alpha_full <- solve(A_lambda, b)
cat("PCG solver correct: ", 
    all.equal(alpha_pcg, alpha_full, tol=1e-10),  "\n")
