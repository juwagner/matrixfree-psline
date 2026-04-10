# ------------------------------------------------------------------------------
# High-level matrix-free operators for additive penalized spline smoothing
#
# This file wraps the low-level tensor‑product and Khatri–Rao operations
# implemented in C++ and provides high-level operators commonly used in
# additive penalized spline regression.
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------
# High-level operations for additive penalized splines
# -----------------------------------------------------------------

library(Rcpp)

sourceCpp("src/base/matrix_free_operations.cpp")
source("src/pspline/pspline_operations.R")

# ------------------------------------------------------------------------------
# Matrix-free multiplication of Φ %*% x, with Φ = sum_s Φ_s
mvp_Phi_terms <- function(PhiT_terms, alpha_terms) {
  n_terms <- length(PhiT_terms)
  v <- rowSums(
    sapply(
      1:n_terms, function(s) 
        mvp_transposed_khatrirao(PhiT_terms[[s]], alpha_terms[[s]])
      )
    )
  return(v)
}

# ------------------------------------------------------------------------------
#Matrix-free multiplication: Φᵀ %*% x, with Φ = sum_s Φ_s
mvp_PhiT_terms <- function(PhiT_terms, x) {
  n_terms <- length(PhiT_terms)
  w <- lapply(1:n_terms, function(s) mvp_khatrirao(PhiT_terms[[s]], x))
  return(w)
}

# ------------------------------------------------------------------------------
# Matrix-free multiplication: (Φᵀ Φ) %*% x, with Φ = sum_s Φ_s
mvp_PhiTPhi_terms <- function(PhiT_terms, alpha_terms) {
  n_terms <- length(PhiT_terms)
  v <- rowSums(
    sapply(
      1:n_terms, 
      function(s) mvp_transposed_khatrirao(PhiT_terms[[s]], alpha_terms[[s]])
      )
    )
  w <- lapply(1:n_terms, function(s) mvp_khatrirao(PhiT_terms[[s]], v))
  return(w)
}

# ------------------------------------------------------------------------------
# Matrix-free multiplication: Λ %*% x, with Λ = blockdiag(Λ_s)
mvp_Lambda_terms <- function(L_terms, alpha_terms) {
  n_terms <- length(L_terms)
  w <- lapply(1:n_terms, function(s) mvp_Lambda(L_terms[[s]], alpha_terms[[s]]))
  return(w)
}

# ------------------------------------------------------------------------------
# Matrix-free multiplication: (Φᵀ Φ + Λ(λ)) %*% x,
# with Φ = sum_s Φ_s and Λ(λ) = blockdiag(λ_s Λ_s)
mvp_lambda_Lambda_terms <- function(L_terms, lambda_vec, alpha_terms) {
  n_terms <- length(L_terms)
  w <- lapply(
    1:n_terms, 
    function(s) lambda_vec[[s]]*mvp_Lambda(L_terms[[s]], alpha_terms[[s]])
    )
  return(w)
}

# -----------------------------------------------
# Matrix-free multiplication: (Φᵀ Φ + Λ(λ)) %*% x
# with Φ = sum_s Φ_s and Λ(λ) = blockdiag(λ_s Λ_s)
mvp_A_lambda_terms <- function(PhiT_terms, L_terms, lambda_vec, alpha_terms) {
  n_terms <- length(L_terms)
  spline <- mvp_PhiTPhi_terms(PhiT_terms, alpha_terms)
  penalty <- mvp_lambda_Lambda_terms(L_terms, lambda_vec, alpha_terms)
  w <- lapply(1:n_terms, function(s) spline[[s]] + penalty[[s]])
  return(w)
}
