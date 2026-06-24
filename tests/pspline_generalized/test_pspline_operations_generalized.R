# ------------------------------------------------------------------------------
# Test suite for high-level P-spline matrix-free operators
# ------------------------------------------------------------------------------

library(microbenchmark)
library(rTensor)

source("src/pspline/pspline_matrices.R")
source("src/pspline_generalized/pspline_operations_generalized.R")

# ------------------------------------------------------------------------------
# Generate test data
# ------------------------------------------------------------------------------

# P-Spline setup
P  <- 3                         # Dimension of tensor-product spline
m  <- c(15, 17, 21)             # interior knots per dimension
q  <- c(3, 2, 3)                # spline degrees
l  <- c(2, 2, 2)                # penalty difference orders
J  <- m + q + 1                 # basis sizes per dimension
n  <- 1000                      # number of data points

# Random input data
X <- matrix(runif(n * P), nrow = n, ncol = P)
alpha <- rnorm(prod(J))
y <- rnorm(n)

# P-Spline bases and penalties
PhiT_list <- lapply(1:P, function(p) {
  build_univarate_bspline_basis_T(X[,p], m[p], q[p])
})

L_list <- lapply(1:P, function(p){
  build_penalty_difference(J[p], l[p])
})

Phi_alpha <- mvp_Phi(PhiT_list, alpha)
W <- exp(Phi_alpha)

# Full matrices for reference
PhiT <- Reduce(rTensor::khatri_rao, PhiT_list)
Phi <- t(PhiT)
PhiT_W_Phi <- PhiT %*% diag(W) %*% Phi

Lambda <- Reduce(`+`, lapply(1:P, function(p) {
  left  <- if(p>1) diag(prod(J[1:(p-1)])) else 1
  right <- if(p<P) diag(prod(J[(p+1):P])) else 1
  kronecker(left, kronecker(L_list[[p]], right))
}))

lambda <- 0.1
A_W_lambda <- PhiT_W_Phi + lambda*Lambda

mvp_ref <- function(A, x) {
  as.vector(A %*% x)
}

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

# mvp_PhiT_W_Phix
mvp_PhiT_W_Phi_mf <- mvp_PhiT_W_Phi(PhiT_list, W, alpha)
mvp_PhiT_W_Phi_ref <- mvp_ref(PhiT_W_Phi, alpha)
cat("MVP with Phi^T %*% W %*% Phi correct:", 
    all.equal(mvp_PhiT_W_Phi_mf, mvp_PhiT_W_Phi_ref, tol=1e-10), "\n")

# mvp_A_W_lambda
mvp_system_mf <- mvp_A_W_lambda(PhiT_list, L_list, W, lambda, alpha)
mvp_system_ref <- mvp_ref(A_W_lambda, alpha)
cat("MVP with weigthed system matrix correct:", 
    all.equal(mvp_system_mf, mvp_system_ref, tol=1e-10), "\n")
