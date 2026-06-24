# ------------------------------------------------------------------------------
# Diagonal-preconditioned Conjugate Gradient (PCG) solver for additive P-splines
# ------------------------------------------------------------------------------

library(Rcpp)

sourceCpp("src/base/matrix_free_operations.cpp")

# ------------------------------------------------------------------------------
# Compute diag(Λ) where Λ = sum_p (I ⊗ L_p ⊗ I)
# L_list = List of penalty matrices L_p (each J_p × J_p)
get_diag_Lambda <- function(L_list) {
  P <- length(L_list)
  J_vec <- vapply(L_list, function(L_p) ncol(L_p), integer(1))
  K <- prod(J_vec)
  
  if (P == 1L) {
    d <- as.vector(diag(L_list[[1L]]))
    return(d)
  }
  
  cp_fwd  <- cumprod(J_vec)
  n_left  <- c(1L, cp_fwd[-P])
  
  cp_bwd  <- cumprod(rev(J_vec))
  n_right <- c(rev(cp_bwd)[-1L], 1L)
  
  diag_pen <- numeric(K)
  for (p in seq_len(P)) {
    dp <- diag(L_list[[p]])
    block <- rep(rep(dp, each = n_right[p]), times = n_left[p])
    diag_pen <- diag_pen + block
  }
  return(as.vector(diag_pen))
}

# ------------------------------------------------------------------------------
# Solve (ΦᵀΦ + Λ(λ)) α = b using diagonal-preconditioned CG
# with Φ = sum_s Φ_s and Λ(λ) = blockdiag(λ_s Λ_s) 
solve_pcg_terms <- function(
    PhiT_terms,
    L_terms,
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
    Aalpha_terms <- mvp_A_lambda_terms(
      PhiT_terms, L_terms, lambda_vec, alpha_terms
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
    Ad_terms <- mvp_A_lambda_terms(PhiT_terms, L_terms, lambda_vec, d_terms)
  
    step_len <- as.numeric(rz / crossprod(unlist(d_terms),unlist(Ad_terms)))
    
    alpha_terms <- lapply(
      1:n_terms, 
      function(s) alpha_terms[[s]] + step_len*d_terms[[s]]
    )
    
    r_terms <- lapply(
      1:n_terms, 
      function(s) r_terms[[s]] - step_len*Ad_terms[[s]]
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

