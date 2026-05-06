# ------------------------------------------------------------------------------
# Diagonal-preconditioned Conjugate Gradient (PCG) solver for P-splines
# ------------------------------------------------------------------------------

library(Rcpp)

sourceCpp("src/base/matrix_free_operations.cpp")
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
    r <- b - mvp_A_w_lambda(PhiT_list, L_list, W, lambda, alpha)
  }
  z <- preconditioner*r
  d <- z
  rz <- as.numeric( crossprod(r,z) )
  
  for(k in 1:it_max){
    Ad <- mvp_A_w_lambda(PhiT_list, L_list, W, lambda, d)
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
