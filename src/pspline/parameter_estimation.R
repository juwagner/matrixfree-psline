# ------------------------------------------------------------------------------
# Estimation of the regularization parameter λ for P-splines
# ------------------------------------------------------------------------------

source("src/pspline/pspline_operations.R")
source("src/pspline/pcg_solver.R")

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
# Estimate df(λ) = trace(S_λ), with S_λ = (A_λ)^{-1} ΦᵀΦ
# Uses Hutchinson trace estimator: trace(S_λ) ≈ 1/M * sum_m  v_mᵀ S_λ v_m,
# where v_m are Rademacher vectors
estimate_trace <- function(
    PhiT_list, L_list, lambda, V_rad, pcg_tol = 10^(-4), pcg_verbose=FALSE
) {
  stopifnot(is.matrix(V_rad))
  K <- nrow(V_rad)
  M <- ncol(V_rad)
  
  trace_terms <- numeric(M)
  for (m in seq_len(M)) {
    v <- V_rad[, m]
    v_tilde <- lambda * mvp_Lambda(L_list=L_list, x=v) 
    v_bar <- solve_pcg(
      PhiT_list=PhiT_list, 
      L_list=L_list, 
      lambda=lambda, 
      b=v_tilde, 
      tol = pcg_tol,
      verbose = pcg_verbose
      )
    trace_terms[m] <- drop(crossprod(v, v_bar))
  }
  df <- K - mean(trace_terms)
  return(as.numeric(df))
}

# ------------------------------------------------------------------------------
# Estimate λ for P-splines using moment-based iteration
# Iteratively updates λ according to:
# λ <- σ²_ε / σ²_α
# where
# σ²_ε = mean((y - Φα)²),
# σ²_α = αᵀ Λ α / df(λ)
estimate_lambda <- function(
    PhiT_list,
    L_list,
    y,
    b,
    lambda_init = 0.1,
    it_max = 10,
    pcg_tol = 1e-3,
    V_rad = NULL,
    M = 5,
    verbose = TRUE
) {
  
  J_vec <- vapply(L_list, ncol, integer(1))
  K <- prod(J_vec)
  
  if (is.null(V_rad)) {
    V_rad <- rademacher_matrix(K, M, seed = 1)
  }
  
  lambda <- lambda_init
  
  alpha <- solve_pcg(
    PhiT_list = PhiT_list,
    L_list = L_list,
    lambda = lambda,
    b = b,
    tol = pcg_tol
  )
  
  for (i in seq_len(it_max)) {
    y_pred <- mvp_Phi(PhiT_list, alpha)
    
    sigma2_eps <- mean((y - y_pred)^2)
    
    df_hat <- estimate_trace(
      PhiT_list = PhiT_list,
      L_list = L_list,
      lambda = lambda,
      V_rad = V_rad,
      pcg_tol = pcg_tol
    )
    
    sigma2_alpha <- drop(crossprod(alpha, mvp_Lambda(L_list, alpha))) / df_hat
    
    lambda_new <- sigma2_eps / sigma2_alpha
    
    if (verbose) {
      cat("Iter ", i, ": lambda = ", lambda_new, " df = ", df_hat, "\n")
    }
    
    if(abs(lambda - lambda_new) <= 0.001) {
      break
    }
    
    lambda <- lambda_new

    alpha <- solve_pcg(
      PhiT_list = PhiT_list,
      L_list = L_list,
      lambda = lambda,
      b = b,
      alpha_init = alpha,
      tol = pcg_tol
    )
    
  }
  return(list(
    lambda = lambda,
    alpha = alpha,
    df = df_hat
  ))
}
