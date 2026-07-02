# ------------------------------------------------------------------------------
# # Test suite: Regularization parameter estimation for generalized P-splines
# ------------------------------------------------------------------------------

library(microbenchmark)

source("src/pspline/pspline_matrices.R")
source("src/pspline_generalized/parameter_estimation_generalized.R")

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

set.seed(42)

#Random input data
X <- matrix(runif(n * P), nrow = n, ncol = P)
y <- as.vector(exp(sin(2*pi*X[,1]*X[,2])*cos(2*pi*X[,3])) + rnorm(n, sd=0.1))

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

# ------------------------------------------------------------------------------
# Estimate alpha for fixed lambda
n_iter <- 5

# With matrix-free algorithm
alpha_mf <- estimate_alpha_generalized(
  n_iter,
  y,
  PhiT_list, 
  L_list,
  lambda,
  alpha_init=alpha,
  pcg_tol=10^(-10),
  pcg_verbose=FALSE
)

# Full computation for comparison
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

cat("Estimation of alpha correct: ", 
    all.equal(alpha_mf, alpha_full, tol=1e-10), "\n")

# ------------------------------------------------------------------------------
# Estimate trace of (ΦᵀWΦ + λΛ)^{-1} λΛ

V_rad <- rademacher_matrix(K, M = 20, seed = 42)
W2 <- as.vector(exp(2 * mvp_Phi(PhiT_list, alpha)))

# With matrix-free algorithm
trace_est <- estimate_trace_generalized(
  PhiT_list = PhiT_list,
  L_list    = L_list,
  W         = W2,
  lambda    = lambda,
  V_rad     = V_rad,
  pcg_tol   = 1e-10
)

# Full computation for comparison
A_W_lambda_true <- PhiT %*% (W2 * Phi) + lambda * Lambda
trace_ref <- sum(diag(solve(A_W_lambda_true) %*% (lambda * Lambda)))

cat("trace estimation <= 5% relative error: ",
    abs(trace_est - trace_ref) / trace_ref <= 0.05, "\n")

# ------------------------------------------------------------------------------
# Estimate lambda for fixed alpha

# With matrix-free algorithm
lambda_est <- estimate_lambda_generalized(
  PhiT_list = PhiT_list,
  L_list    = L_list,
  alpha     = alpha,
  lambda    = lambda,
  V_rad     = V_rad,
  pcg_tol   = 1e-10
)

# Full computation for comparison
#W2 <- as.vector(exp(2 * mvp_Phi(PhiT_list, alpha)))
#A_W_lambda_true <- PhiT %*% (W2 * Phi) + lambda * Lambda
#trace_ref <- sum(diag(solve(A_W_lambda_true) %*% (lambda * Lambda)))
W1 <- as.vector(exp(mvp_Phi(PhiT_list, alpha)))
sigma_eps_ref   <- mean((y - W1)^2)
sigma_alpha_ref <- drop(crossprod(alpha, mvp_Lambda(L_list, alpha))) / (K - trace_ref)
lambda_ref      <- sigma_eps_ref / sigma_alpha_ref

cat("lambda estimation <= 5% relative error: ",
    abs(lambda_est - lambda_ref) / lambda_ref <= 0.05, "\n")

# ------------------------------------------------------------------------------
# Estimate alpha and lambda in parallel

result <- fit_pspline_generalized(
  n_iter        = 3,
  n_iter_alpha  = 3,
  y             = y,
  PhiT_list     = PhiT_list,
  L_list        = L_list,
  lambda        = lambda,
  V_rad         = V_rad,
  pcg_tol       = 1e-4
)

#print(result)

alpha_hat <- result$alpha
y_hat <- exp(mvp_Phi(PhiT_list = PhiT_list, x = alpha_hat))

cat("RMSE of y and y_hat: ", sqrt(mean((y-y_hat)^2)), "\n")
cat("RRMSE of y and y_hat: ", sqrt(mean((y-y_hat)^2)) / mean(y), "\n")
