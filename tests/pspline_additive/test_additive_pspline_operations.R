# ------------------------------------------------------------------
# Testing functions
# ------------------------------------------------------------------

library(rTensor)
library(Matrix)

source("src/pspline/pspline_matrices.R")
source("src/additive_pspline/additive_pspline_operations.R")

# ------------------------------------------------------------
# Random test data

# Setup
n_summands <- 2

P  <- list(3, 2)
m  <- list( c(11, 7, 12), c(8, 10) )
q <- list( c(3, 2, 3), c(3, 2) )
l <- list( c(2, 2, 2), c(2, 2) )
J <- lapply(1:n_summands, function(s) m[[s]]+q[[s]]+1)
K <- lapply(1:n_summands, function(s) prod(J[[s]]))
n  <- 5000

# Data
X_terms <- lapply(1:n_summands, function(s) matrix(runif(n * P[[s]]), nrow = n, ncol = P[[s]]))
alpha_terms <- lapply(1:n_summands, function(s) rnorm(prod(J[[s]])))
y <- rnorm(n)

# P-Spline
PhiT_terms <- lapply(1:n_summands, function(s) lapply(1:P[[s]], function(p) {
  build_univarate_bspline_basis_T(X_terms[[s]][,p], m[[s]][p], q[[s]][p])
}))

L_terms <- lapply(1:n_summands, function(s) lapply(1:P[[s]], function(p){
  build_penalty_difference(J[[s]][p], l[[s]][p])
}))

# Full matrices for P-Splines
PhiT_full <- lapply(1:n_summands, function(s) Reduce(rTensor::khatri_rao, PhiT_terms[[s]]))
PhiT_full <- do.call(rbind, PhiT_full)
Phi <- t(PhiT_full)
PhiTPhi_full <- PhiT_full %*% Phi

Lambda_list <- lapply(1:n_summands, function(s){
  Reduce(`+`, lapply(1:P[[s]], function(p) {
    left  <- if(p>1) diag(prod(J[[s]][1:(p-1)])) else 1
    right <- if(p<P[[s]]) diag(prod(J[[s]][(p+1):P[[s]]])) else 1
    kronecker(left, kronecker(L_terms[[s]][[p]], right))
  }))
})

lambda_vec <- list(0.1, 0.2)

Lambda_list_weighted <- lapply(1:n_summands, function(s) lambda_vec[[s]] * Lambda_list[[s]])
lambda_Lambda_full <- bdiag(Lambda_list_weighted)

A_lambda_full <- PhiTPhi_full + lambda_Lambda_full

alpha <- unlist(alpha_terms)

# ------------------------------------------------------------
# Reference function in R
# ------------------------------------------------------------

mvp_ref <- function(A, x) {
  as.vector(A %*% x)
}

# ------------------------------------------------------------
# Correctness test
# ------------------------------------------------------------

cat("=== Correctness Test ===\n")

# Phi
mvp_Phi_mf <- mvp_Phi_terms(PhiT_terms, alpha_terms)
mvp_Phi_ref <- mvp_ref(Phi, alpha)
cat("MVP with Phi correct:", all.equal(mvp_Phi_mf, mvp_Phi_ref, tol=1e-10), "\n")

# PhiT
mvp_PhiT_mf <- unlist(mvp_PhiT_terms(PhiT_terms, y))
mvp_PhiT_ref <- mvp_ref(PhiT_full, y)
cat("MVP with Phi^T correct:", all.equal(mvp_PhiT_mf, mvp_PhiT_ref, tol=1e-10), "\n")

# mvp_PhiTPhialpha
mvp_PhiTPhi_mf <- unlist(mvp_PhiTPhi_terms(PhiT_terms, alpha_terms))
mvp_PhiTPhi_ref <- mvp_ref(PhiTPhi_full, alpha)
cat("MVP with Phi^T %*% Phi correct:", all.equal(mvp_PhiTPhi_mf, mvp_PhiTPhi_ref, tol=1e-10), "\n")

# mvp_Lambda
mvp_Lambda_mf <- unlist(mvp_lambda_Lambda_terms(L_terms, lambda_vec, alpha_terms))
mvp_Lambda_ref <- mvp_ref(lambda_Lambda_full, alpha)
cat("MVP with (weighted) Lambda correct:", all.equal(mvp_Lambda_mf, mvp_Lambda_ref, tol=1e-10), "\n")

# mvp_A_lambda
mvp_A_lambda_mf <- unlist(mvp_A_lambda_terms(PhiT_terms, L_terms, lambda_vec, alpha_terms))
mvp_A_lambda_ref <- mvp_ref(A_lambda_full, alpha)
cat("MVP with A_lambda correct:", all.equal(mvp_A_lambda_mf, mvp_A_lambda_ref, tol=1e-10), "\n")
