# ------------------------------------------------------------------------------
# High-level matrix-free operators for generalized penalized spline smoothing
#
# This file wraps the low-level tensor‑product and Khatri–Rao operations
# implemented in C++ and provides high-level operators commonly used in
# generalized penalized spline regression.
# ------------------------------------------------------------------------------

library(Rcpp)

sourceCpp("src/base/matrix_free_operations.cpp")
source("src/pspline/pspline_operations.R")

# ------------------------------------------------------------------------------
# Matrix-free multiplication: (Φᵀ W Φ) %*% x, where W is a vector representing
# the diagonal of a weight matrix.
mvp_PhiT_W_Phi <- function(PhiT_list, W, x) {
  Phix <- mvp_Phi(PhiT_list, x)
  return(mvp_PhiT(PhiT_list, W*Phix))
}

# ------------------------------------------------------------------------------
# Matrix-free multiplication: (Φᵀ W Φ + λ Λ) %*% x
mvp_A_w_lambda <- function(PhiT_list, L_list, lambda, x) {
  return(mvp_PhiT_W_Phi(PhiT_list, x) + lambda * mvp_Lambda(L_list, x))
}
