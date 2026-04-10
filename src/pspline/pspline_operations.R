# ------------------------------------------------------------------------------
# High-level matrix-free operators for penalized spline smoothing
#
# This file wraps the low-level tensor‑product and Khatri–Rao operations
# implemented in C++ and provides high-level operators commonly used in
# penalized spline regression.
# ------------------------------------------------------------------------------

library(Rcpp)

sourceCpp("src/base/matrix_free_operations.cpp")

# ------------------------------------------------------------------------------
# Matrix-free multiplication: Φ %*% x
mvp_Phi <- function(PhiT_list, x) {
  return(mvp_transposed_khatrirao(PhiT_list, x))
}

# ------------------------------------------------------------------------------
#Matrix-free multiplication: Φᵀ %*% x
mvp_PhiT <- function(PhiT_list, x) {
  return(mvp_khatrirao(PhiT_list, x))
}

# ------------------------------------------------------------------------------
# Matrix-free multiplication: (Φᵀ Φ) %*% x
mvp_PhiTPhi <- function(PhiT_list, x) {
  return(mvp_gram_khatrirao(PhiT_list, x))
}

# ------------------------------------------------------------------------------
# Matrix-free multiplication: Λ %*% x
mvp_Lambda <- function(L_list, x) {
  
  P <- length(L_list)
  J_vec <- sapply(1:P, function(p) dim(L_list[[p]])[1])
  n_left <- c(1, sapply(1:(P-1), function(p) prod(J_vec[1:p]) ))
  n_right <- c(rev( sapply(1:(P-1), function(p) prod(rev(J_vec)[1:p]) ) ), 1)
  if(P==1) {
    n_left <- n_right <- 1
  }
  Lambda_x <- rowSums(
    sapply(
      1:P, function(p) mvp_normalfactor(L_list[[p]], n_left[p], n_right[p], x)
      )
    )
  return(as.vector(Lambda_x))
  
}

# ------------------------------------------------------------------------------
# Matrix-free multiplication: (Φᵀ Φ + λ Λ) %*% x
mvp_A_lambda <- function(PhiT_list, L_list, lambda, x) {
  return(mvp_gram_khatrirao(PhiT_list, x) + lambda * mvp_Lambda(L_list, x))
}
