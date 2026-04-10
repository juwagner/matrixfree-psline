# -----------------------------------------------------------------
# High-level operations for additive penalized splines
# -----------------------------------------------------------------

library(Rcpp)

sourceCpp("src/matrix_free_operations.cpp")
source("R/pspline/pspline_operations.R")

# -----------------------------------------------
mvp_Phi_terms <- function(PhiT_terms, alpha_terms) {
  n_terms <- length(PhiT_terms)
  v <- rowSums(sapply(1:n_terms, function(s) mvp_transposed_khatrirao(PhiT_terms[[s]], alpha_terms[[s]])))
  return(v)
}

# -----------------------------------------------
mvp_PhiT_terms <- function(PhiT_terms, y) {
  n_terms <- length(PhiT_terms)
  w <- lapply(1:n_terms, function(s) mvp_khatrirao(PhiT_terms[[s]], y))
  return(w)
}

# -----------------------------------------------
mvp_PhiTPhi_terms <- function(PhiT_terms, alpha_terms) {
  n_terms <- length(PhiT_terms)
  v <- rowSums(sapply(1:n_terms, function(s) mvp_transposed_khatrirao(PhiT_terms[[s]], alpha_terms[[s]])))
  w <- lapply(1:n_terms, function(s) mvp_khatrirao(PhiT_terms[[s]], v))
  return(w)
}

# -----------------------------------------------
mvp_Lambda_terms <- function(L_terms, alpha_terms) {
  n_terms <- length(L_terms)
  w <- lapply(1:n_terms, function(s) mvp_Lambda(L_terms[[s]], alpha_terms[[s]]))
  return(w)
}

# -----------------------------------------------
mvp_lambda_Lambda_terms <- function(L_terms, lambda_vec, alpha_terms) {
  n_terms <- length(L_terms)
  w <- lapply(1:n_terms, function(s) lambda_vec[[s]]*mvp_Lambda(L_terms[[s]], alpha_terms[[s]]))
  return(w)
}

# -----------------------------------------------
mvp_A_lambda_terms <- function(PhiT_terms, L_terms, lambda_vec, alpha_terms) {
  n_terms <- length(L_terms)
  spline <- mvp_PhiTPhi_terms(PhiT_terms, alpha_terms)
  penalty <- mvp_lambda_Lambda_terms(L_terms, lambda_vec, alpha_terms)
  w <- lapply(1:n_terms, function(s) spline[[s]] + penalty[[s]])
  return(w)
}
