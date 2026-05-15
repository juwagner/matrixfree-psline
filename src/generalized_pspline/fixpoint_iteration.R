# ------------------------------------------------------------------------------
# Diagonal-preconditioned Conjugate Gradient (PCG) solver for P-splines
# ------------------------------------------------------------------------------

source("src/pspline/pcg_solver.R")
source("src/generalized_pspline/generalized_pspline_operations.R")

# ------------------------------------------------------------------------------
# ToDo: compute diag(Φᵀ W Φ) as preconditioner (e.g. using Hutchinson)

# ------------------------------------------------------------------------------
# Solve (Φᵀ W Φ + λΛ) α = b using diagonal-preconditioned CG
solve_generalized_pcg = function(
    PhiT_list, 
    L_list,
    W,
    lambda, 
    b,
    alpha_init=NULL, 
    it_max=length(b),
    tol=10^(-4), 
    verbose=FALSE
){
  
  P <- length(L_list)
  J_vec <- sapply(1:P, function(p) dim(L_list[[p]])[2] )
  K <- prod(J_vec)
  
  diag_PhiTPhi <- get_diag_PhiTPhi(PhiT_list = PhiT_list)
  diag_Lambda <- get_diag_Lambda(L_list = L_list)
  preconditioner <- 1 / (diag_PhiTPhi + lambda*diag_Lambda)
  
  norm_b <- sqrt(drop(crossprod(b)))
  
  if(is.null(alpha_init)){
    r <- b
    alpha <- rep(0,K)
  } else{
    alpha <- alpha_init
    r <- b - mvp_A_W_lambda(PhiT_list, L_list, W, lambda, alpha)
  }
  z <- preconditioner*r
  d <- z
  rz <- as.numeric( crossprod(r,z) )
  
  for(k in 1:it_max){
    Ad <- mvp_A_W_lambda(PhiT_list, L_list, W, lambda, d)
    step_len <- as.numeric(rz / crossprod(d,Ad))
    alpha <- alpha + step_len*d
    r <- r - step_len*Ad
    z <- preconditioner*r
    rz_old <- rz
    rz <- as.numeric(crossprod(r,z))
    
    relres <- sqrt(drop(crossprod(r))) / norm_b
    if(verbose == TRUE) {
      cat("PCG iteration: ", k, " relres: ", relres , "\n" )
    }
      
    if(relres < tol){
      break
    }
    
    beta <- rz / rz_old
    d <- z + beta*d
    
  }
  
  return(as.vector(alpha))
  
}

# ------------------------------------------------------------------------------
# Run fixpoint iteration to solve for alpha in generalized pspline
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
    
    a <- lambda * mvp_Lambda(L_list, alpha)
    z <- mvp_PhiT(PhiT_list , W1*(y-W1)) - a
    v <- solve_generalized_pcg(
      PhiT_list, L_list, W2, lambda, z, alpha, verbose=pcg_verbose, tol=pcg_tol
    )
    
    alpha <- alpha + v
  }
  
  return(as.vector(alpha))

}


