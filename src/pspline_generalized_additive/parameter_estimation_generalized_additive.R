# ------------------------------------------------------------------------------
# Estimation of the regularization parameters λ_s for generalized additive P-splines
# ------------------------------------------------------------------------------

source("src/pspline_additive/pspline_operations_additive.R")
source("src/pspline_additive/parameter_estimation_additive.R")
source("src/pspline_generalized_additive/pspline_operations_generalized_additive.R")
source("src/pspline_generalized_additive/pcg_solver_generalized_additive.R")


# ------------------------------------------------------------------------------
# Iteration to estimate α_terms (with fixed λ_s) in generalized additive p-spline model
estimate_alpha_generalized_additive = function(
    n_iter,
    PhiT_terms, 
    L_terms,
    lambda_vec, 
    alpha_init=NULL,
    pcg_tol=10^(-4),
    pcg_verbose=FALSE
){
  n_terms <- length(PhiT_terms)
  if(is.null(alpha_init)){
    J_vec_terms <- lapply(L_terms, function(Ls) vapply(Ls, ncol, 1L))
    K_terms <- vapply(
      L_terms, function(Ls) prod(vapply(Ls, ncol, numeric(1))), numeric(1)
    )
    alpha_terms <- lapply(seq_len(n_terms), function(s) rep(0, K_terms[s]))
  } else{
    alpha_terms <- alpha_init
  }
  
  for(i in 1:n_iter){
    cat("---------- Start fixpoint iteration: ", i, "\n")
    Phi_alpha <- mvp_Phi_terms(PhiT_terms, alpha_terms)
    W1 <- as.vector(exp(Phi_alpha))
    W2 <- as.vector(exp(2*Phi_alpha))
    
    rhs <- mvp_PhiT_terms(PhiT_terms , W1*(y-W1)) - mvp_lambda_Lambda_terms(L_terms, lambda_vec, alpha_terms)
    v <- solve_pcg_generalized_terms(
      PhiT_terms = PhiT_terms,
      L_terms = L_terms,
      W = W2,
      lambda_vew = lambda_vec,
      b_terms = rhs,
      alpha_init = alpha, 
      verbose=pcg_verbose,
      tol=pcg_tol
    )
    
    alpha_new <- lapply(1:n_terms, function(p) alpha_terms[[p]] + v[[p]])
    rel <- mean((unlist(alpha_terms) - unlist(alpha_new))^2)
    cat("Relative change of alpha: ", rel, "\n")
    alpha_terms <- alpha_new
  }
  
  return(as.vector(alpha_terms))
  
}



