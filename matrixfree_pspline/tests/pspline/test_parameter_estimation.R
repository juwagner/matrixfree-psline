# ------------------------------------------------------------------------------
# Test suite: Regularization parameter estimation for P-splines
# ------------------------------------------------------------------------------

library(rTensor)

source("src/pspline/pspline_matrices.R")
source("src/pspline/parameter_estimation.R")

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

# Random input data
X <- matrix(runif(n * P), nrow = n, ncol = P)
y <- as.vector(sin(2*pi*X[,1]*X[,2]) + 0.3*cos(2*pi*X[,2]) + rnorm(n, sd=0.1))

# P-Spline bases and penalties
PhiT_list <- lapply(1:P, function(p) {
  build_univarate_bspline_basis_T(X[,p], m[p], q[p])
})

L_list <- lapply(1:P, function(p){
  build_penalty_difference(J_vec[p], l[p])
})

# Full matrices for reference
PhiT <- Reduce(rTensor::khatri_rao, PhiT_list)
Phi <- t(PhiT)
PhiTPhi <- PhiT %*% Phi

Lambda <- Reduce(`+`, lapply(1:P, function(p) {
  left  <- if(p>1) diag(prod(J_vec[1:(p-1)])) else 1
  right <- if(p<P) diag(prod(J_vec[(p+1):P])) else 1
  kronecker(left, kronecker(L_list[[p]], right))
}))

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Hutchinson DF estimation
lambda <- 0.01
A_lambda <- PhiTPhi + lambda*Lambda
S_lambda <- solve(A_lambda) %*% PhiTPhi

df_ref <- sum(diag(S_lambda))

V_rad <- rademacher_matrix(K=K, M=100, seed=42)
df_est <- estimate_df(PhiT_list, L_list, lambda, V_rad)

cat("DF estimation <= 5% relative error: ", 
    (abs(df_est - df_ref) / df_ref) <= 0.05,"\n")

# ------------------------------------------------------------------------------
# lambda estimation
V_rad <- rademacher_matrix(K=K, M=10, seed=42)

estimation <- estimate_lambda(
  PhiT_list = PhiT_list,
  L_list = L_list,
  y = y,
  b = mvp_PhiT(PhiT_list, y),
  lambda_init = 1,
  it_max = 10,
  verbose = TRUE
)

print(estimation)


