# ------------------------------------------------------------------------------
# Estimation of α and λ for generalized P-splines
# ------------------------------------------------------------------------------

source("src/utils/rademacher.R")
source("src/pspline_generalized/pspline_operations_generalized.R")
source("src/pspline_generalized/pcg_solver_generalized.R")

# ------------------------------------------------------------------------------
# Iteration to estimate α (with fixed λ) in generalized p-spline model
estimate_alpha_generalized = function(
    n_iter,
    y,
    PhiT_list,
    L_list,
    lambda,
    alpha_init=NULL,
    pcg_tol=10^(-4),
    pcg_verbose=FALSE
){
  if(is.null(alpha_init)){
    P <- length(L_list)
    J_vec <- sapply(1:P, function(p) dim(L_list[[p]])[2] )
    K <- prod(J_vec)
    alpha <- rep(0,K)
  } else{
    alpha <- alpha_init
  }
  
  for(i in 1:n_iter){
    cat("Iteration for alpha: ", i, "/", n_iter, "\n")
    Phi_alpha <- mvp_Phi(PhiT_list, alpha)
    W1 <- as.vector(exp(Phi_alpha))
    W2 <- as.vector(exp(2*Phi_alpha))
    
    rhs <- mvp_PhiT(PhiT_list , W1*(y-W1)) - lambda * mvp_Lambda(L_list, alpha)
    v <- solve_pcg_generalized(
      PhiT_list, L_list, W2, lambda, rhs, alpha, verbose=pcg_verbose, tol=pcg_tol
    )
    
    alpha_new <- alpha + v
    rel <- mean((alpha-alpha_new)^2)
    cat("Relative change of alpha: ", rel, "\n")
    alpha <- alpha_new
  }
  
  cat("Solved for alpha \n")
  return(as.vector(alpha))
  
}

# ------------------------------------------------------------------------------
# Estimate trace of (ΦᵀWΦ + λΛ)^{-1} λΛ using Hutchinson.
estimate_trace_generalized = function(
    PhiT_list, 
    L_list,
    W,
    lambda,
    V_rad,
    pcg_tol = 10^(-4), 
    pcg_verbose=FALSE
){
  stopifnot(is.matrix(V_rad))
  K <- nrow(V_rad)
  M <- ncol(V_rad)
  
  trace_terms <- numeric(M)
  for (m in seq_len(M)) {
    v <- V_rad[, m]
    v_tilde <- lambda * as.vector(mvp_Lambda(L_list=L_list, x=v))
    v_bar <- solve_pcg_generalized(
      PhiT_list=PhiT_list, 
      L_list=L_list,
      W=W,
      lambda=lambda,
      b=v_tilde, 
      tol = pcg_tol,
      verbose = pcg_verbose
    )
    trace_terms[m] <- drop(crossprod(v, v_bar))
  }
  
  cat("Estimated trace \n")
  return(as.numeric(mean(trace_terms)))
  
}

# ------------------------------------------------------------------------------
# Iteration to estimate λ (for fixed α) in generalized p-spline model
estimate_lambda_generalized = function(
    PhiT_list, 
    L_list,
    alpha,
    lambda=0.1,
    V_rad,
    pcg_tol = 10^(-4), 
    pcg_verbose=FALSE
  ){
  
  K <- length(alpha)
  Phi_alpha <- mvp_Phi(PhiT_list, alpha)
  W1 <- as.vector(exp(Phi_alpha))
  W2 <- as.vector(exp(2*Phi_alpha))
  
  sigma_eps <- mean((y-W1)^2)
  
  trace_est <- estimate_trace_generalized(PhiT_list, L_list, W2, lambda, V_rad)
  sigma_alpha <- crossprod(alpha, mvp_Lambda(L_list, alpha)) / (K - trace_est)
  
  lambda <- as.numeric(sigma_eps / sigma_alpha)
  
  cat("Estimated lambda \n")
  return(lambda)
}


# ------------------------------------------------------------------------------
# Iteration to estimate α and λ in parallel in generalized p-spline model
fit_pspline_generalized = function(
    n_iter=2,
    n_iter_alpha=2,
    y,
    PhiT_list,
    L_list,
    lambda,
    V_rad,
    pcg_tol=10^(-2)
    ){

  alpha <- NULL

  for (i in 1:n_iter) {
    cat("--- Outer iteration: ", i ,"/", n_iter, "\n")
    alpha <- estimate_alpha_generalized(
      n_iter=n_iter_alpha,
      y=y,
      PhiT_list=PhiT_list,
      L_list=L_list,
      lambda=lambda,
      alpha_init=alpha,
      pcg_tol=pcg_tol
    )

    if(i != n_iter){
      lambda <- estimate_lambda_generalized(
        PhiT_list=PhiT_list, 
        L_list=L_list,
        alpha=alpha,
        lambda=lambda,
        V_rad=V_rad,
        pcg_tol=pcg_tol
      )
      cat("Current lambda: ", lambda, "\n")
    }
  }
  
  return(list(
    lambda = lambda,
    alpha = alpha
  ))
  
}
