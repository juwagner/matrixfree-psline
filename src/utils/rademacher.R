# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Generate Rademacher random matrix
rademacher_matrix <- function(K, M, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  V_rad <- matrix(
    sample(c(-1L, 1L), size = K * M, replace = TRUE), nrow = K, ncol = M
  )
  storage.mode(V_rad) <- "double"
  return(V_rad)
}

# ------------------------------------------------------------------------------
# Generate Rademacher random matrix per term
rademacher_matrix_terms <- function(K_terms, M, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  V <- lapply(
    K_terms, 
    function(K) matrix(
      sample(c(-1L, 1L), size = K*M, replace = TRUE), nrow = K, ncol = M
    )
  )
  return(V)
}