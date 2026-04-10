# ------------------------------------------------------------------------------
# Spline construction + difference penalty matrices
#
# These functions are used to construct:
# - Univariate B-spline bases using truncated power functions
# - Corresponding transposed basis matrices
# - Difference penalty matrices (Dᵀ D) of arbitrary order
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Build univariate B-spline basis (truncated power basis formulation)
build_univarate_bspline_basis <- function(x, m, q, Omega=c(min(x),max(x))) {
  
  if(m < 1) stop("Number of interior knots 'm' must be >= 1")
  if(q < 0) stop("Spline degree 'q' must be >= 0")
  
  h <- (Omega[2]-Omega[1]) / (m+1)
  knots <- seq(Omega[1]-q*h, Omega[2]+q*h, by=h)
  
  dx <- outer(x, knots, "-")
  H <- pmax(dx, 0)^q
  
  K <- length(knots)
  D <- diff(diag(K), diff=q+1) / (gamma(q+1) * h^q)
  
  Phi <- (-1)^(q+1) * tcrossprod(H, D)
  
  Phi[abs(Phi) < 1e-10] <- 0
  
  return(Phi)
}

# ------------------------------------------------------------------------------
# Build transposed univariate B-spline basis (convenience wrapper)
build_univarate_bspline_basis_T<- function(x, m, q, Omega=c(min(x),max(x))) {
  return(t(build_univarate_bspline_basis(x, m, q, Omega=c(min(x),max(x)))))
}

# ------------------------------------------------------------------------------
# Build difference penalty matrix (Dᵀ D)
build_penalty_difference <- function(J, l = 2) {
  
  if(J <= l) stop("Dimension J must be larger than difference order l")
  
  diff_matrix <- diff(diag(J), differences = l)
  L <- crossprod(diff_matrix)
  return(L)
  
}