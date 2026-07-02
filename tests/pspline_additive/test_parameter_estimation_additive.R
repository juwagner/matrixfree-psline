# ------------------------------------------------------------------------------
# Test suite: Regularization parameter estimation for additive P-splines
# ------------------------------------------------------------------------------

rm(list=ls())

library(rTensor)
library(Matrix)

source("src/pspline/pspline_matrices.R")
source("src/pspline_additive/parameter_estimation_additive.R")

# ------------------------------------------------------------
# Random test data
# ------------------------------------------------------------

source("tests/pspline_additive/generate_test_data_additive.R")

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

V_rad_terms <- rademacher_matrix_terms(K_terms=K_terms, M=10, seed=42)

# ------------------------------------------------------------------------------
# df estimation
S_lambda <- solve(A_lambda_full) %*% PhiTPhi_full
df_ref <- sum(diag(S_lambda))

df_est <- estimate_df_terms(PhiT_terms, L_terms, lambda_vec, V_rad_terms)

cat("df estimation <= 5% relative error: ", 
    (abs(df_est - df_ref) / df_ref) <= 0.05,"\n")

# ------------------------------------------------------------------------------
# lambda estimation
estimation <- estimate_lambda_terms(
  PhiT_terms = PhiT_terms,
  L_terms = L_terms,
  y = y,
  b_terms = mvp_PhiT_terms(PhiT_terms, y),
  V_rad_terms = V_rad_terms,
  verbose = TRUE
)

print(estimation)
