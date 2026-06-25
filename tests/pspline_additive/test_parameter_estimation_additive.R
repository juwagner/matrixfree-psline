# ------------------------------------------------------------------------------
# Test suite: Regularization parameter estimation for additive P-splines
# ------------------------------------------------------------------------------

library(rTensor)
library(Matrix)

source("src/pspline/pspline_matrices.R")
source("src/pspline_additive/parameter_estimation_additive.R")

# ------------------------------------------------------------
# Random test data
# ------------------------------------------------------------

# Additive P-Spline setup
n_terms <- 2

P  <- list(3, 2)                    # Dimension of spline per term
m  <- list(c(11, 7, 12), c(8, 10))  # interior knots per term per dimension
q <- list(c(3, 2, 3), c(3, 2))      # spline degrees per term
l <- list(c(2, 2, 2), c(2, 2))      # penalty difference orders per term
J_vec <- lapply(
  1:n_terms,
  function(s) m[[s]]+q[[s]]+1       # basis sizes per dimension
)
K_terms <- lapply(J_vec, function(s) prod(s))
n  <- 5000                          # number of data points


# Random input data
X_terms <- lapply(
  1:n_terms, 
  function(s) matrix(runif(n * P[[s]]), nrow = n, ncol = P[[s]])
)
alpha_terms <- lapply(
  1:n_terms, 
  function(s) rnorm(prod(J_vec[[s]]))
)
y <- sin(2*pi*X_terms[[1]][,1]*X_terms[[1]][,2])*cos(2*pi*X_terms[[1]][,3]) + 
  0.3*cos(2*pi*X_terms[[2]][,1]*X_terms[[2]][,2]) + rnorm(n, sd=0.1)

# P-Spline bases and penalties
PhiT_terms <- lapply(
  1:n_terms, 
  function(s) lapply(
    1:P[[s]], function(p) 
      build_univarate_bspline_basis_T(X_terms[[s]][,p], m[[s]][p], q[[s]][p])
  )
)

L_terms <- lapply(
  1:n_terms, 
  function(s) lapply(1:P[[s]], function(p)
    build_penalty_difference(J_vec[[s]][p], l[[s]][p])
  )
)

# Full matrices for reference
PhiT_full <- lapply(
  1:n_terms, function(s) Reduce(rTensor::khatri_rao, PhiT_terms[[s]])
)
PhiT_full <- do.call(rbind, PhiT_full)
Phi <- t(PhiT_full)
PhiTPhi_full <- PhiT_full %*% Phi

Lambda_list <- lapply(1:n_terms, function(s){
  Reduce(`+`, lapply(1:P[[s]], function(p) {
    left  <- if(p>1) diag(prod(J_vec[[s]][1:(p-1)])) else 1
    right <- if(p<P[[s]]) diag(prod(J_vec[[s]][(p+1):P[[s]]])) else 1
    kronecker(left, kronecker(L_terms[[s]][[p]], right))
  }))
})

lambda_vec <- c(0.1, 0.2)

Lambda_list_weighted <- lapply(1:n_terms, function(s) lambda_vec[[s]] * Lambda_list[[s]])
lambda_Lambda_full <- bdiag(Lambda_list_weighted)

A_lambda_full <- PhiTPhi_full + lambda_Lambda_full

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

V_rad_terms <- rademacher_matrix_terms(K_terms=K_terms, M=10, seed=42)

# ------------------------------------------------------------------------------
# trace estimation

# already tested

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
