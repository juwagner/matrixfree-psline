# ------------------------------------------------------------------------------
# Test suite for PCG solver + diagonal matrix-free operators for P-splines
# ------------------------------------------------------------------------------

library(microbenchmark)

source("src/pspline/pspline_matrices.R")
source("src/pspline_generalized/pcg_solver_generalized.R")

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
y <- as.vector(rnorm(n))

# P-Spline bases and penalties
PhiT_list <- lapply(1:P, function(p) {
  build_univarate_bspline_basis_T(X[,p], m[p], q[p])
})

L_list <- lapply(1:P, function(p){
  build_penalty_difference(J_vec[p], l[p])
})

alpha <- rnorm(prod(J_vec))
W2 <- exp(2*mvp_Phi(PhiT_list, alpha))

# Full matrices for reference
PhiT <- Reduce(rTensor::khatri_rao, PhiT_list)
Phi <- t(PhiT)
PhiT_W_Phi <- PhiT %*% diag(W2) %*% Phi

Lambda <- Reduce(`+`, lapply(1:P, function(p) {
  left  <- if(p>1) diag(prod(J_vec[1:(p-1)])) else 1
  right <- if(p<P) diag(prod(J_vec[(p+1):P])) else 1
  kronecker(left, kronecker(L_list[[p]], right))
}))

lambda <- 0.1
A_W_lambda <- PhiT_W_Phi + lambda*Lambda

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

# PCG solver
b <- as.vector(PhiT %*% y)
alpha_pcg  <- solve_pcg_generalized(PhiT_list, L_list, W2, lambda, b, tol = 1e-12, verbose=TRUE)
alpha_full <- solve(A_W_lambda, b)
cat("PCG solver correct: ", 
    all.equal(alpha_pcg, alpha_full, tol=1e-10),  "\n")
