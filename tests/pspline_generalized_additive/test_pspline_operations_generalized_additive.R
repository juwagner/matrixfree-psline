# ------------------------------------------------------------------------------
# Test suite for high-level additive P-spline matrix-free operators
# ------------------------------------------------------------------------------

library(rTensor)
library(Matrix)

source("src/pspline/pspline_matrices.R")
source("src/pspline_generalized_additive/pspline_operations_generalized_additive.R")

# ------------------------------------------------------------------------------
# Generate test data
# ------------------------------------------------------------------------------

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
n  <- 5000                          # number of data points

# Random input data
X_terms <- lapply(
  1:n_terms, 
  function(s) matrix(runif(n * P[[s]]), nrow = n, ncol = P[[s]])
)
alpha_terms <- lapply(1:n_terms, function(s) rnorm(prod(J_vec[[s]])))
y <- rnorm(n)

# P-Spline bases and penalties
PhiT_terms <- lapply(
  1:n_terms, 
  function(s) 
    lapply(1:P[[s]], function(p)
      build_univarate_bspline_basis_T(X_terms[[s]][,p], m[[s]][p], q[[s]][p])
  )
)

L_terms <- lapply(
  1:n_terms, 
  function(s) lapply(1:P[[s]], function(p)
    build_penalty_difference(J_vec[[s]][p], l[[s]][p])
  )
)

Phi_terms_alpha <- mvp_Phi_terms(PhiT_terms, alpha_terms)
W <- exp(Phi_terms_alpha)

# Full matrices for reference
PhiT_full <- lapply(
  1:n_terms, 
  function(s) Reduce(rTensor::khatri_rao, PhiT_terms[[s]])
)
PhiT_full <- do.call(rbind, PhiT_full)
Phi <- t(PhiT_full)
PhiT_W_Phi_full <- PhiT_full %*% diag(W) %*% Phi

Lambda_list <- lapply(1:n_terms, function(s){
  Reduce(`+`, lapply(1:P[[s]], function(p) {
    left  <- if(p>1) diag(prod(J_vec[[s]][1:(p-1)])) else 1
    right <- if(p<P[[s]]) diag(prod(J_vec[[s]][(p+1):P[[s]]])) else 1
    kronecker(left, kronecker(L_terms[[s]][[p]], right))
  }))
})

lambda_vec <- list(0.1, 0.2)

Lambda_list_weighted <- lapply(
  1:n_terms, 
  function(s) lambda_vec[[s]] * Lambda_list[[s]]
)
lambda_Lambda_full <- bdiag(Lambda_list_weighted)

A_W_lambda_full <- PhiT_W_Phi_full + lambda_Lambda_full

alpha <- unlist(alpha_terms)

mvp_ref <- function(A, x) {
  as.vector(A %*% x)
}

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

# mvp_PhiT_W_Phi_terms
mvp_PhiT_W_Phi_mf <- unlist(mvp_PhiT_W_Phi_terms(PhiT_terms, W, alpha_terms))
mvp_PhiT_W_Phi_ref <- mvp_ref(PhiT_W_Phi_full, alpha)
cat("MVP with Phi^T %*% W %*% Phi correct:", 
    all.equal(mvp_PhiT_W_Phi_mf, mvp_PhiT_W_Phi_ref, tol=1e-10), "\n")

# mvp_A_W_lambda_terms
mvp_system_mf <- unlist(mvp_A_W_lambda_terms(PhiT_terms, L_terms, lambda_vec, W, alpha_terms))
mvp_system_ref <- mvp_ref(A_W_lambda_full, alpha)
cat("MVP with weigthed system matrix correct:", 
    all.equal(mvp_system_mf, mvp_system_ref, tol=1e-10), "\n")