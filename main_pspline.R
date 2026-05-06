# ------------------------------------------------------------------------------
# P-spline models for the LUCAS dataset
# ------------------------------------------------------------------------------

rm(list=ls())

library(Rcpp)
library(tictoc)

sourceCpp("src/base/matrix_free_operations.cpp")
source("src/pspline/pspline_matrices.R")
source("src/pspline/pspline_operations.R")
source("src/pspline/pcg_solver.R")
source("src/pspline/parameter_estimation.R")

# ------------------------------------------------------------------------------
# load data
source("src/utils/load_lucas_data.R")

X_input <- X           # U
main <- "y = s(x)"     # s(u)

# ------------------------------------------------------------------------------
# P-spline setup

P     <- ncol(X_input)        # dimension of tensor product
m     <- rep(36, P)           # interior knots per dimension
q     <- rep(3,  P)           # spline degree per dimension
l     <- rep(2,  P)           # difference penalty order
J_vec <- m + q + 1            # basis sizes
K     <- prod(J_vec)          # total number of coefficients


PhiT_list <- lapply(
  1:P, function(p) build_univarate_bspline_basis_T(X_input[,p], m[p], q[p])
)

L_list <- lapply(1:P, function(p) build_penalty_difference(J=J_vec[p], l=l[p]))

b <- mvp_PhiT(PhiT_list = PhiT_list, x = y)

# ------------------------------------------------------------------------------
# Solve for α using a fixed λ (single PCG run)

lambda <- 0.01184

tic("single PCG run")

alpha <- solve_pcg(
  PhiT_list = PhiT_list,
  L_list = L_list,
  lambda = lambda,
  b = b, 
  verbose=TRUE
)

toc()

y_hat <- mvp_Phi(PhiT_list, alpha)

# ------------------------------------------------------------------------------
# Estimate λ and α simultaneously

V_rad <- rademacher_matrix(K=K, M=3, seed=42) 

tic("total time for lambda estimation")

estimation <- estimate_lambda(
  PhiT_list = PhiT_list,
  L_list = L_list,
  y = y,
  b = b,
  lambda_init = 0.05,
  it_max = 10,
  V_rad = V_rad,
)

toc()

alpha  <- estimation$alpha
lambda <- estimation$lambda
df     <- estimation$df

y_hat <- mvp_Phi(PhiT_list, alpha)

# ------------------------------------------------------------------------------
# Validation metrics

res <- y-y_hat
RSS <- sum(res^2)
V_rad <- rademacher_matrix(K=K, M=3, seed=42)
df <- estimate_trace(PhiT_list, L_list, lambda, V_rad)
AIC <- 2*n*log(RSS) + 2*df

cat(
  "P-Spline Model Validation | RSS:", RSS, 
  "DF:", df, 
  "AIC:", AIC, 
  "Min fitted:", min(y_hat), 
  "\n"
)
 
# ------------------------------------------------------------------------------
# Diagnostic plots

qqnorm(res, main=main)
qqline(res)
plot(y_hat, res, xlab="Fitted Values", ylab="Residuals", main=main)

