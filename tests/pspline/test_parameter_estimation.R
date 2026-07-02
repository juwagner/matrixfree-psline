# ------------------------------------------------------------------------------
# Test suite: Regularization parameter estimation for P-splines
# ------------------------------------------------------------------------------

rm(list=ls())

library(rTensor)

source("src/pspline/pspline_matrices.R")
source("src/pspline/parameter_estimation.R")

# ------------------------------------------------------------------------------
# Generate test data
# ------------------------------------------------------------------------------

source("tests/pspline/generate_test_data.R")

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Hutchinson DF estimation
lambda <- 0.01
A_lambda <- PhiTPhi + lambda*Lambda
S_lambda <- solve(A_lambda) %*% PhiTPhi

df_ref <- sum(diag(S_lambda))

V_rad <- rademacher_matrix(K=K, M=100, seed=42)
df_est <- estimate_trace(PhiT_list, L_list, lambda, V_rad)

cat("DF estimation <= 5% relative error: ", 
    (abs(df_est - df_ref) / df_ref) <= 0.05,"\n")

# ------------------------------------------------------------------------------
# lambda estimation
V_rad <- rademacher_matrix(K=K, M=10, seed=42)

estimation <- estimate_lambda(
  PhiT_list = PhiT_list,
  L_list = L_list,
  y = y,
  b = mvp_PhiT(PhiT_list, y),
  lambda_init = 1,
  it_max = 10,
  verbose = TRUE
)

print(estimation)


