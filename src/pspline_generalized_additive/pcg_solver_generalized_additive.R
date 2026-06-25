# ------------------------------------------------------------------------------
# Diagonal-preconditioned Conjugate Gradient (PCG) solver for generalized additive P-splines
# ------------------------------------------------------------------------------

library(Rcpp)

sourceCpp("src/base/matrix_free_operations.cpp")
source("src/pspline_additive/pcg_solver_additive.R")
source("src/pspline_generalized_additive/pspline_operations_generalized_additive.R")

# ------------------------------------------------------------------------------
# Solve (Φᵀ W Φ + Λ(λ)) α = b using diagonal-preconditioned CG
# with Φ = sum_s Φ_s and Λ(λ) = blockdiag(λ_s Λ_s) 
solve_pcg_generalized_terms <- function(
    PhiT_terms,
    L_terms,
    W,
    lambda_vec,
    b_terms,
    alpha_init = NULL,
    it_max = NULL,
    tol = 1e-4,
    verbose = FALSE
) {
  
  if(is.null(it_max)) {
    it_max <- sum(
      vapply(
        L_terms, function(Ls) prod(vapply(Ls, ncol, numeric(1))), numeric(1)
      )
    )
  }
  
  n_terms <- length(PhiT_terms)
  J_vec_terms <- lapply(L_terms, function(Ls) vapply(Ls, ncol, 1L))
  K_terms <- vapply(
    L_terms, function(Ls) prod(vapply(Ls, ncol, numeric(1))), numeric(1)
  )
  diag_PhiTPhi_terms <- lapply(PhiT_terms, diag_gram_khatrirao)
  diag_Lambda_terms  <- lapply(L_terms, get_diag_Lambda)
  preconditioner_terms <- lapply(seq_len(n_terms), function(s) {
    1 / (diag_PhiTPhi_terms[[s]] + lambda_vec[s] * diag_Lambda_terms[[s]])
  })

  norm_b <- sqrt(sum(unlist(b_terms)^2))
  if (is.null(alpha_init)) {
    alpha_terms <- lapply(seq_len(n_terms), function(s) rep(0, K_terms[s]))
    r_terms <- b_terms
  } else {
    alpha_terms <- alpha_init
    A_W_alpha_terms <- mvp_A_W_lambda_terms(
      PhiT_terms, L_terms, lambda_vec, W, alpha_terms
    )
    r_terms <- lapply(1:n_terms, function(s) b_terms[[s]] - Aalpha_terms[[s]])
  }
  z_terms <- lapply(
    1:n_terms, 
    function(s) preconditioner_terms[[s]] * r_terms[[s]]
  )
  d_terms <- z_terms
  rz <- as.numeric(crossprod(unlist(r_terms), unlist(z_terms)))
  
  for (k in seq_len(it_max)) {
    A_W_d_terms <- mvp_A_W_lambda_terms(PhiT_terms, L_terms, lambda_vec, W, d_terms)
  
    step_len <- as.numeric(rz / crossprod(unlist(d_terms),unlist(A_W_d_terms)))
    
    alpha_terms <- lapply(
      1:n_terms, 
      function(s) alpha_terms[[s]] + step_len*d_terms[[s]]
    )
    
    r_terms <- lapply(
      1:n_terms, 
      function(s) r_terms[[s]] - step_len*A_W_d_terms[[s]]
    )
    
    relres <- sqrt(sum(unlist(r_terms)^2)) / norm_b
    if (verbose) {
      cat("PCG iter:", k, "relres:", relres, "\n")
    }
    if (relres < tol) break
    
    z_terms <- lapply(
      1:n_terms, 
      function(s) preconditioner_terms[[s]]*r_terms[[s]]
    )
    
    rz_new <- as.numeric(crossprod(unlist(r_terms), unlist(z_terms)))
    beta <- rz_new / rz
    rz <- rz_new
    
    d_terms <- lapply(1:n_terms, function(s) z_terms[[s]] + beta*d_terms[[s]])
  }
  return(alpha_terms)
}

