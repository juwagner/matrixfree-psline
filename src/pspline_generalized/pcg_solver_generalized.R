# ------------------------------------------------------------------------------
# Diagonal-preconditioned Conjugate Gradient (PCG) solver for generalized P-splines
# ------------------------------------------------------------------------------

source("src/pspline_generalized/pspline_operations_generalized.R")

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
# Compute diag(ΦᵀΦ) matrix-free
# PhiT_list List of transposed marginal B-spline bases Φ_pᵀ.
get_diag_PhiTPhi <- function(PhiT_list) {
  diag_gram_khatrirao(A_list = PhiT_list)
}

# ------------------------------------------------------------------------------
# ToDo: compute diag(Φᵀ W Φ) as preconditioner

# ------------------------------------------------------------------------------
# Solve (Φᵀ W Φ + λΛ) α = b using diagonal-preconditioned CG
solve_pcg_generalized = function(
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
