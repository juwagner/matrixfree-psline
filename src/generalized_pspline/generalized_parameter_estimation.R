# ------------------------------------------------------------------------------
# Estimation of α and λ for generalized P-splines
# ------------------------------------------------------------------------------

source("src/generalized_pspline/generalized_pspline_operations.R")
source("src/generalized_pspline/generalized_pcg_solver.R")

# ------------------------------------------------------------------------------
# Fixpoint iteration to solve for α in generalized pspline
fixpoint_iteration_alpha = function(
    n_iter,
    PhiT_list, 
    L_list,
    lambda, 
    b,
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
    cat("---------- Start fixpoint iteration: ", i,"\n" )
    Phi_alpha <- mvp_Phi(PhiT_list, alpha)
    W1 <- exp(Phi_alpha)
    W2 <- exp(2*Phi_alpha)
    
    score <- mvp_PhiT(PhiT_list , W1*(y-W1)) - lambda * mvp_Lambda(L_list, alpha)
    v <- solve_generalized_pcg(
      PhiT_list, L_list, W2, lambda, score, alpha, verbose=pcg_verbose, tol=pcg_tol
    )
    
    alpha <- alpha + v
  }
  
  return(as.vector(alpha))
  
}
