# ------------------------------------------------------------------------------
# Test suite for matrix-free tensor and Khatri–Rao operations
#
# This script validates correctness and benchmarks performance of the Rcpp
# implementations in `matrix_free_operations.cpp`.
# ------------------------------------------------------------------------------

library(Rcpp)
library(microbenchmark)
library(rTensor)

sourceCpp("src/base/matrix_free_operations.cpp")

# ------------------------------------------------------------------------------
# Generate test data
# ------------------------------------------------------------------------------

# Normalfactor test data
J     <- 20     # dimension of the middle matrix A
left  <- 30     # left identity dimension
right <- 40     # right identity dimension

J     <- 20
left  <- 30
right <- 40

A_inner_nf <- matrix(rnorm(J*J), J, J)
A_nf <- kronecker(diag(left), kronecker(A_inner_nf, diag(right)))
x_nf <- rnorm(left * J * right)

# Khatri–Rao test data
P  <- 3                      # number of matrices in Khatri–Rao product
m  <- rep(20, P)             # row sizes of A₁, A₂, A₃
n  <- 10000                  # number of columns in each A_p

A_list_kr <- lapply(m, function(mm) matrix(rnorm(mm * n), mm, n))

# Explicit dense Khatri–Rao matrix for reference
A_kr <- Reduce(rTensor::khatri_rao, A_list_kr)
A_transposed_kr <- t(A_kr)
A_gram_kr <- A_kr %*% A_transposed_kr

x_kr <- rnorm(n)
y_kr <- rnorm(prod(m))

mvp_ref <- function(A, x) {
  as.vector(A %*% x)
}

# ------------------------------------------------------------------------------
# Correctness tests
# ------------------------------------------------------------------------------

cat("======= Correctness Test =======\n")

# Normalfactor
res_nf_mf <- mvp_normalfactor(A_inner_nf, left, right, x_nf)
res_nf_ref <- mvp_ref(A_nf, x_nf)
cat("Normalfactor correct:", all.equal(res_nf_mf, res_nf_ref, tol=1e-10), "\n")

# Khatri-Rao matrix
res_kr_mf <- mvp_khatrirao(A_list_kr, x_kr)
res_kr_ref <- mvp_ref(A_kr, x_kr)
cat("Khatri-Rao matrix product correct:", 
    all.equal(res_kr_mf, res_kr_ref, tol=1e-10), "\n")

# Transposed Khatri-Rao matrix
res_tkr_mf <- mvp_transposed_khatrirao(A_list_kr, y_kr)
res_tkr_ref <- mvp_ref(A_transposed_kr, y_kr)
cat("Transposed Khatri-Rao matrix product correct:", 
    all.equal(res_tkr_mf, res_tkr_ref, tol=1e-10), "\n")

# Gram matrix multiplication
res_gram_mf <- mvp_gram_khatrirao(A_list_kr, y_kr)
res_gram_ref <- mvp_ref(A_gram_kr, y_kr)
cat("Gram Khatri-Rao matrix product correct:", 
    all.equal(res_gram_mf, res_gram_ref, tol=1e-10), "\n")

# Diagonal of Gram matrix
res_diag_mf <- diag_gram_khatrirao(A_list_kr)
res_diag_ref <- diag(A_gram_kr)
cat("Diagonal of gram Khatri-Rao matrix correct:", 
    all.equal(res_diag_mf, res_diag_ref, tol=1e-10), "\n")

# ------------------------------------------------------------------------------
# Performance benchmark
# ------------------------------------------------------------------------------

cat("======= Performance Benchmark =======\n")

res <- microbenchmark(
  
  normalfactor_mf = mvp_normalfactor(A_inner_nf, left, right, x_nf),
  normalfactor_ref = mvp_ref(A_nf, x_nf),
  
  khatrirao_mf = mvp_khatrirao(A_list_kr, x_kr),
  khatrirao_ref = mvp_ref(A_kr, x_kr),
  
  transposed_khatrirao_mf = mvp_transposed_khatrirao(A_list_kr, y_kr),
  transposed_khatrirao_ref = mvp_ref(A_transposed_kr, y_kr),
  
  gram_khatrirao_mf = mvp_gram_khatrirao(A_list_kr, y_kr),
  gram_khatrirao_ref = mvp_ref(A_gram_kr, y_kr),
  
  diag_gram_khatrirao_mf = diag_gram_khatrirao(A_list_kr),
  diag_gram_khatrirao_ref = diag(A_gram_kr),
  
  times = 10,
  unit = "ms"
)

print(res)