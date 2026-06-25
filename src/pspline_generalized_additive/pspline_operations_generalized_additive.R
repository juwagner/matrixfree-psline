# ------------------------------------------------------------------------------
# High-level matrix-free operators for generalized additive penalized spline smoothing
#
# This file wraps the low-level tensor‑product and Khatri–Rao operations
# implemented in C++ and provides high-level operators commonly used in
# generalized additive penalized spline regression.
# ------------------------------------------------------------------------------

library(Rcpp)

sourceCpp("src/base/matrix_free_operations.cpp")
source("src/pspline/pspline_operations.R")
source("src/pspline_additive/pspline_operations_additive.R")

# ------------------------------------------------------------------------------
# Matrix-free multiplication: (Φᵀ W Φ) %*% x, with Φ = sum_s Φ_s
mvp_PhiT_W_Phi_terms <- function(PhiT_terms, W, alpha_terms) {
  x <- mvp_Phi_terms(PhiT_terms=PhiT_terms, alpha_terms=alpha_terms)
  return(mvp_PhiT_terms(PhiT_terms=PhiT_terms, x=W*x))
}

# -----------------------------------------------
# Matrix-free multiplication: (Φᵀ W Φ + Λ(λ)) %*% x
# with Φ = sum_s Φ_s and Λ(λ) = blockdiag(λ_s Λ_s)
mvp_A_W_lambda_terms <- function(PhiT_terms, L_terms, lambda_vec, W, alpha_terms) {
  n_terms <- length(L_terms)
  spline <- mvp_PhiT_W_Phi_terms(PhiT_terms, W, alpha_terms)
  penalty <- mvp_lambda_Lambda_terms(L_terms, lambda_vec, alpha_terms)
  w <- lapply(1:n_terms, function(s) spline[[s]] + penalty[[s]])
  return(w)
}
