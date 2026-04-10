# -----------------------------------------------------------------
# Estimation of the regularization parameter
# -----------------------------------------------------------------

source("R/pspline/pcg_solver.R")
source("R/pspline_additive/additive_pspline_operations.R")
source("R/pspline_additive/additive_pcg_solver.R")

# -----------------------------------------------
rademacher_terms <- function(K_terms, M, seed = NUL) {
  if (!is.null(seed)) {
    set.seed(seed)
  } 
  n_terms <- length(K_terms)
  V <- lapply(1:M, function(m) lapply(1:n_terms, function(s) sample(c(-1,1), K_terms[s], replace=TRUE)))
  return(V)
}

# -----------------------------------------------
estimate_trace_terms <- function(
    PhiT_terms, L_terms, lambda_vec, V_rad_terms, pcg_tol = 1e-4, pcg_verbose=FALSE
) {
  n_terms <- length(PhiT_terms)
  K_terms <- vapply(V_rad_terms[[1]], function(Ms) length(Ms), numeric(1))
  M <- length(V_rad_terms)
  V1 <- lapply(1:M, function(m) mvp_lambda_Lambda_terms(L_terms, lambda_vec, V_rad_terms[[m]]) ) 
  V2 <- lapply(1:M, function(m) lapply(1:n_terms, function(s) solve_pcg_single(PhiT_terms[[s]], L_terms[[s]], lambda_vec[s], V1[[m]][[s]]) ) )
  trace <- sapply(1:n_terms, function(s) K_terms[s] - mean( sapply(1:M, function(m) crossprod( V_rad_terms[[m]][[s]],V2[[m]][[s]] ) ) ) )
  return(trace)
}

# -----------------------------------------------
estimate_lambda_terms <- function(
    PhiT_terms,
    L_terms,
    y,
    b_terms,
    lambda_vec_init = NULL,
    it_max = 10,
    pcg_tol = 1e-4,
    M = 5,
    seed = NULL,
    V_rad_terms = NULL,
    verbose = TRUE
) {
  
  n_terms <- length(PhiT_terms)
  
  K_terms <- vapply(L_terms, function(Ls) prod(vapply(Ls, ncol, numeric(1))), numeric(1))
  
  if (is.null(V_rad_terms)) {
    if (!is.null(seed)) set.seed(seed)
    V_rad_terms <- rademacher_matrix_terms(K_terms, M)
  }
  
  if (is.null(lambda_vec_init)) {
    lambda_vec <- rep(0.1, n_terms)
  } else {
    lambda_vec <- lambda_vec_init
  }
  
  alpha_terms <- solve_pcg_terms(
    PhiT_terms, L_terms, lambda_vec, b_terms,
    tol = pcg_tol
  )
  
  for (i in seq_len(it_max)) {
    
    f_terms <- lapply(seq_len(n_terms), function(s) {
      mvp_Phi(PhiT_terms[[s]], alpha_terms[[s]])
    })
    
    y_pred <- Reduce("+", f_terms)
    
    sigma2_eps <- mean((y - y_pred)^2)
    
    trace_hat <- estimate_trace_terms(
      PhiT_terms = PhiT_terms,
      L_terms    = L_terms,
      lambda_vec = lambda_vec,
      V_rad_terms    = V_rad_terms,
      pcg_tol    = pcg_tol
    )
    
    sigma2_alpha_terms <- vapply(
      seq_len(n_terms),
      function(s) {
        drop(crossprod(alpha_terms[[s]],
                       mvp_Lambda(L_terms[[s]], alpha_terms[[s]]))) / trace_hat[s]
      },
      numeric(1)
    )
    
    lambda_new <- sigma2_eps / sigma2_alpha_terms
    
    if (verbose) {
      cat("Iter", i, ": lambda =", paste(round(lambda_new, 6), collapse = " , "), "\n")
    }
    
    if (max(abs(lambda_new - lambda_vec)) < 0.001)
      break
    
    lambda_vec <- lambda_new
    
    alpha_terms <- solve_pcg_terms(
      PhiT_terms, L_terms, lambda_vec, b_terms, alpha_init=alpha_terms,
      tol = pcg_tol
    )
  }
  
  return(list(
    lambda_vec = lambda_vec,
    alpha_terms = alpha_terms,
  ))
}

# -----------------------------------------------
# estimate_df_terms <- function(PhiT_terms, L_terms, lambda_vec, V_rad_terms, pcg_tol = 1e-4, pcg_verbose=FALSE) {
#   n_terms <- length(PhiT_terms)
#   K_terms <- vapply(V_rad_terms, function(Ms) nrow(Ms), numeric(1))
#   M <- ncol(V_rad_terms[[1]])
#   trace_terms <- numeric(M)
#   for (m in seq_len(M)) {
#     v_terms <- lapply(seq_len(n_terms), function(s) V_rad_terms[[s]][, m])
#     w_terms <- mvp_lambda_Lambda_terms(L_terms, lambda_vec, v_terms)
#     u_terms <- solve_pcg_terms(
#       PhiT_terms, L_terms, lambda_vec, w_terms, tol = pcg_tol, verbose = pcg_verbose
#     )
#     trace_terms[m] <- sum(sum(unlist(v_terms) * unlist(u_terms)))
#     
#   }
#   df_total <- sum(K_terms) - mean(trace_terms)
#   return(as.numeric(df_total))
# }