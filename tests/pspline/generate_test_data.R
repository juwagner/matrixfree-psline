# ------------------------------------------------------------------------------
# Generate data for P-spline tests
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
alpha <- rnorm(K)

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

lambda <- 0.1
A_lambda <- PhiTPhi + lambda*Lambda