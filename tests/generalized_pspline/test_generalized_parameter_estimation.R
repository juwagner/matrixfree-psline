# ------------------------------------------------------------------------------
# Test suite for PCG solver + diagonal matrix-free operators for P-splines
# ------------------------------------------------------------------------------

library(microbenchmark)

source("src/pspline/pspline_matrices.R")
source("src/generalized_pspline/generalized_parameter_estimation.R")

# ------------------------------------------------------------------------------
# Generate test data
# ------------------------------------------------------------------------------

# P-Spline setup
P  <- 3                         # Dimension of tensor-product spline
m  <- c(15, 11, 12)             # interior knots per dimension
q  <- c(3, 2, 3)                # spline degrees
l  <- c(2, 2, 2)                # penalty difference orders
J_vec  <- m + q + 1             # basis sizes per dimension
K <- prod(J_vec)                # overall basis size
n  <- 5000                      # number of data points

#Random input data
X <- matrix(runif(n * P), nrow = n, ncol = P)
y <- as.vector(sin(2*pi*X[,1]*X[,2])*cos(2*pi*X[,3]) + rnorm(n, sd=0.1))

sin(2*pi*X[,1]*X[,2])*cos(2*pi*X[,3])

# P-Spline bases and penalties
PhiT_list <- lapply(1:P, function(p) {
  build_univarate_bspline_basis_T(X[,p], m[p], q[p])
})

L_list <- lapply(1:P, function(p){
  build_penalty_difference(J_vec[p], l[p])
})

alpha <- rnorm(prod(J_vec))
Phi_alpha <-mvp_Phi(PhiT_list, alpha)
W1 <- exp(Phi_alpha)
W2 <- exp(2*Phi_alpha)

# Full matrices for reference
PhiT <- Reduce(rTensor::khatri_rao, PhiT_list)
Phi <- t(PhiT)

Lambda <- Reduce(`+`, lapply(1:P, function(p) {
  left  <- if(p>1) diag(prod(J_vec[1:(p-1)])) else 1
  right <- if(p<P) diag(prod(J_vec[(p+1):P])) else 1
  kronecker(left, kronecker(L_list[[p]], right))
}))

lambda <- 0.01

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

# Single iteration
n_iter <- 3

alpha_fixpoint <- fixpoint_w_alpha(
  n_iter,
  PhiT_list, 
  L_list,
  lambda,
  alpha_init=alpha,
  pcg_tol=10^(-6),
  pcg_verbose=FALSE
)

alpha_full <- alpha
for (i in 1:n_iter){
  PhiT <- Reduce(rTensor::khatri_rao, PhiT_list)
  Phi <- t(PhiT)
  Phi_alpha <- Phi %*% alpha_full
  W1 <- as.vector(exp(Phi_alpha))
  W2 <- as.vector(exp(2*Phi_alpha))
  PhiT_W_Phi <- PhiT %*% diag(W2) %*% Phi
  A_W_lambda <- PhiT_W_Phi + lambda*Lambda
  s <- PhiT%*%(diag(W1)%*%(y- W1)) - lambda*Lambda%*%alpha_full
  alpha_full <- as.vector(alpha_full + solve(A_W_lambda)%*%s)
}

cat("Fixpoint iteration correct: ", 
    all.equal(alpha_fixpoint, alpha_full, tol=1e-10),  "\n")
