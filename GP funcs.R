# Functions to sample from and predict values of a Gaussian process.

# Squared distance function
sq_dist <- function(x1, x2) (x1 - x2)^2

# Covariance function - squared exponential
K <- function(x1, x2, l, sigma_f) {
  sigma_f^2 * exp(-outer(x1, x2, sq_dist) / (2 * l^2))
}

# Function to sample from the GP
samp_GP <- function(x, l, sigma_f, n, sigma_n) {
  Sigma <- K(x, x, l, sigma_f)
  L <- try(chol(Sigma), T)
  if (inherits(L, "try-error")) return(numeric(n))
  L %*% rnorm(n)
}

# Function for Moore-Penrose pseudo-inverse via singular value decomposition
pinv <- function(mat) {
  svd_mat <- svd(mat)
  pinv_svs <- ifelse(svd_mat$d == 0, 0, 1 / svd_mat$d)
  svd_mat$v %*% diag(pinv_svs) %*% t(svd_mat$u)
}

# Function to predict from GP
pred_GP <- function(x, x_star, l, sigma_f, sigma_n, n, y) {
  k_star <- K(x, x_star, l, sigma_f)
  k_inv <- pinv(K(x, x, l, sigma_f) + sigma_n^2 * diag(n))
  f_mean = t(k_star) %*% k_inv %*% y
  f_var = K(x_star, x_star, l, sigma_f) - t(k_star) %*% k_inv %*% k_star
  list(f_mean, f_var)
}