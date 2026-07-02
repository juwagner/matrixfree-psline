# ------------------------------------------------------------------------------
# Additive P-spline models for the LUCAS dataset
# ------------------------------------------------------------------------------

rm(list=ls())

library(Rcpp)
library(tictoc)

sourceCpp("src/base/matrix_free_operations.cpp")
source("src/pspline/pspline_matrices.R")
source("src/pspline_additive/pspline_operations_additive.R")
source("src/pspline_additive/pcg_solver_additive.R")
source("src/pspline_additive/parameter_estimation_additive.R")

# ------------------------------------------------------------------------------
# load data

source("src/utils/load_lucas_data.R")

X_terms <- list(X, U)

# ------------------------------------------------------------------------------
# Additive P-spline setup

n_terms <- length(X_terms)
P <- sapply(1:n_terms, function(s) dim(X_terms[[s]])[2])
m <- lapply(1:n_terms, function(s) rep(36,P[s]))
q <- lapply(1:n_terms, function(s) rep(3,P[s]))
l <- lapply(1:n_terms, function(s) rep(2,P[s]))
J <- lapply(1:n_terms, function(s) m[[s]]+q[[s]]+1)
K <- sapply(1:n_terms, function(s) prod(J[[s]]))

PhiT_terms <- lapply(1:n_terms, function(s) lapply(1:P[[s]], function(p) {
  build_univarate_bspline_basis_T(X_terms[[s]][,p], m[[s]][p], q[[s]][p])
}))

L_terms <- lapply(1:n_terms, function(s) lapply(1:P[[s]], function(p){
  build_penalty_difference(J[[s]][p], l[[s]][p])
}))

b_terms = mvp_PhiT_terms(PhiT_terms = PhiT_terms, x = y)

# ------------------------------------------------------------------------------
# Solve for α using a fixed λ (single PCG run)

lambda_vec <- c(0.01, 0.45)

tic("single PCG run")

alpha_terms <- solve_pcg_terms(
  PhiT_terms = PhiT_terms,
  L_terms = L_terms,
  lambda_vec = lambda_vec,
  b_terms = b_terms,
  verbose = TRUE,
)

toc()

y_hat <- mvp_Phi_terms(PhiT_terms = PhiT_terms, alpha_terms = alpha_terms)

# ------------------------------------------------------------------------------
# Estimate λ_vec and α_terms simultaneously

V_rad_terms <- rademacher_matrix_terms(K_terms=K, M=3, seed=42)

tic("total time for lambda estimation")

estimation <- estimate_lambda_terms(
  PhiT_terms = PhiT_terms,
  L_terms = L_terms,
  y = y,
  b_terms = b_terms,
  lambda_vec_init = c(0.1, 0.1),
  it_max = 5,
  V_rad_terms = V_rad_terms,
  verbose = TRUE
)

toc()

alpha_terms <- estimation$alpha_terms
lambda_vec <- estimation$lambda_vec

y_hat <- mvp_Phi_terms(PhiT_terms = PhiT_terms, alpha_terms = alpha_terms)

# ------------------------------------------------------------------------------
# Validation metrics

res <- y-y_hat
RSS <- sum(res^2)
V_rad_terms <- rademacher_matrix_terms(K_terms=K, M=3, seed=42)
df <- estimate_df_terms(PhiT_terms, L_terms, lambda_vec, V_rad_terms)
AIC <- 2*n*log(RSS) + 2*df

cat(
  "P-Spline Model Validation | RSS:", RSS, 
  "DF:", df, 
  "AIC:", AIC, 
  "Min fitted:", min(y_hat), 
  "\n"
)

