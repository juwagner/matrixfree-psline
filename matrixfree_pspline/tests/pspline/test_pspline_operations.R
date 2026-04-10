# ------------------------------------------------------------------------------
# Test suite for high-level P-spline matrix-free operators
# ------------------------------------------------------------------------------

library(microbenchmark)
library(rTensor)

source("src/pspline/pspline_matrices.R")
source("src/pspline/pspline_operations.R")

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

# Full matrices for reference
PhiT <- Reduce(rTensor::khatri_rao, PhiT_list)
Phi <- t(PhiT)
PhiTPhi <- PhiT %*% Phi

Lambda <- Reduce(`+`, lapply(1:P, function(p) {
  left  <- if(p>1) diag(prod(J[1:(p-1)])) else 1
  right <- if(p<P) diag(prod(J[(p+1):P])) else 1
  kronecker(left, kronecker(L_list[[p]], right))
}))

lambda <- 0.1
A_lambda <- PhiTPhi + lambda*Lambda

mvp_ref <- function(A, x) {
  as.vector(A %*% x)
}

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

# Phi
mvp_Phi_mf <- mvp_Phi(PhiT_list, alpha)
mvp_Phi_ref <- mvp_ref(Phi, alpha)
cat("MVP with Phi correct:", 
    all.equal(mvp_Phi_mf, mvp_Phi_ref, tol=1e-10), "\n")

# PhiT
mvp_PhiT_mf <- mvp_PhiT(PhiT_list, y)
mvp_PhiT_ref <- mvp_ref(PhiT, y)
cat("MVP with Phi^T correct:", 
    all.equal(mvp_PhiT_mf, mvp_PhiT_ref, tol=1e-10), "\n")

# mvp_PhiTPhix
mvp_PhiTPhi_mf <- mvp_PhiTPhi(PhiT_list, alpha)
mvp_PhiTPhi_ref <- mvp_ref(PhiTPhi, alpha)
cat("MVP with Phi^T %*% Phi correct:", 
    all.equal(mvp_PhiTPhi_mf, mvp_PhiTPhi_ref, tol=1e-10), "\n")

# mvp_Lambda
mvp_Lambda_mf <- mvp_Lambda(L_list, alpha)
mvp_Lambda_ref <- mvp_ref(Lambda, alpha)
cat("MVP with Lambda correct:", 
    all.equal(mvp_Lambda_mf, mvp_Lambda_ref, tol=1e-10), "\n")

# mvp_A_lambda
mvp_system_mf <- mvp_A_lambda(PhiT_list, L_list, lambda, alpha)
mvp_system_ref <- mvp_ref(A_lambda, alpha)
cat("MVP with system matrix correct:", 
    all.equal(mvp_system_mf, mvp_system_ref, tol=1e-10), "\n")

# ------------------------------------------------------------------------------
# Performance test
# ------------------------------------------------------------------------------

cat("======= Performance Benchmark =======\n")

res <- microbenchmark(
  mvp_Phi_mf = mvp_Phi(PhiT_list, alpha),
  mvp_Phi_ref = mvp_ref(Phi, alpha),
  
  mvp_PhiT_mf = mvp_PhiT(PhiT_list, y),
  mvp_PhiT_ref = mvp_ref(PhiT, y),
  
  mvp_PhiTPhi_mf = mvp_PhiTPhi(PhiT_list, alpha),
  mvp_PhiTPhi_ref = mvp_ref(PhiTPhi, alpha),
  
  mvp_Lambda_mf = mvp_Lambda(L_list, alpha),
  mvp_Lambda_ref = mvp_ref(Lambda, alpha),
  
  mvp_A_lambda_mf = mvp_A_lambda(PhiT_list, L_list, lambda, alpha),
  mvp_A_lambda_ref = mvp_ref(A_lambda, alpha),
  
  times = 10,
  unit = "ms"
)

print(res)