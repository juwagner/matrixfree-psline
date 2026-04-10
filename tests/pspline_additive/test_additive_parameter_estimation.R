# ------------------------------------------------------------------
# Testing functions
# ------------------------------------------------------------------

library(rTensor)
library(Matrix)

source("src/pspline/pspline_matrices.R")
source("src/additive_pspline/additive_parameter_estimation.R")

# ------------------------------------------------------------
# Random test data
# ------------------------------------------------------------

# Setup
n_terms <- 2

P  <- list(3, 2)
m  <- list( c(11, 7, 12), c(8, 10) )
q <- list( c(3, 2, 3), c(3, 2) )
l <- list( c(2, 2, 2), c(2, 2) )
J <- lapply(1:n_terms, function(s) m[[s]]+q[[s]]+1)
K <- lapply(1:n_terms, function(s) prod(J[[s]]))
n  <- 5000

# Data
X_terms <- lapply(1:n_terms, function(s) matrix(runif(n * P[[s]]), nrow = n, ncol = P[[s]]))
alpha_terms <- lapply(1:n_terms, function(s) rnorm(prod(J[[s]])))
y <- sin(2*pi*X_terms[[1]][,1]*X_terms[[1]][,2])*cos(2*pi*X_terms[[1]][,3]) + 0.3*cos(2*pi*X_terms[[2]][,1]*X_terms[[2]][,2]) + rnorm(n, sd=0.1)

# P-Spline
PhiT_terms <- lapply(1:n_terms, function(s) lapply(1:P[[s]], function(p) {
  build_univarate_bspline_basis_T(X_terms[[s]][,p], m[[s]][p], q[[s]][p])
}))

L_terms <- lapply(1:n_terms, function(s) lapply(1:P[[s]], function(p){
  build_penalty_difference(J[[s]][p], l[[s]][p])
}))

# Full matrices for P-Splines
PhiT_full <- lapply(1:n_terms, function(s) Reduce(rTensor::khatri_rao, PhiT_terms[[s]]))
PhiT_full <- do.call(rbind, PhiT_full)
Phi <- t(PhiT_full)
PhiTPhi_full <- PhiT_full %*% Phi

Lambda_list <- lapply(1:n_terms, function(s){
  Reduce(`+`, lapply(1:P[[s]], function(p) {
    left  <- if(p>1) diag(prod(J[[s]][1:(p-1)])) else 1
    right <- if(p<P[[s]]) diag(prod(J[[s]][(p+1):P[[s]]])) else 1
    kronecker(left, kronecker(L_terms[[s]][[p]], right))
  }))
})

lambda_vec <- c(0.1, 0.2)

Lambda_list_weighted <- lapply(1:n_terms, function(s) lambda_vec[[s]] * Lambda_list[[s]])
lambda_Lambda_full <- bdiag(Lambda_list_weighted)

A_lambda_full <- PhiTPhi_full + lambda_Lambda_full

# ------------------------------------------------------------
# Correctness test
# ------------------------------------------------------------

# -----------------------------------------------
# Trace / DF estimation

S_lambda <- solve(A_lambda_full) %*% PhiTPhi_full
df_ref <- sum(diag(S_lambda))

V_rad_terms <- rademacher_matrix_terms(K_terms=K, M=10, seed=42)

df_est <- estimate_df_terms(PhiT_terms, L_terms, lambda_vec, V_rad_terms)

cat("Trace estimation <= 5% relative error: ", (abs(df_est - df_ref) / df_ref) <= 0.05,"\n")

# -----------------------------------------------
# lambda estimation

V_rad_terms <- rademacher_matrix_terms(K_terms=K, M=10, seed=42)

estimation <- estimate_lambda_terms(
  PhiT_terms = PhiT_terms,
  L_terms = L_terms,
  y = y,
  b_terms = mvp_PhiT_terms(PhiT_terms, y),
  V_rad_terms = V_rad_terms,
  verbose = TRUE
)

