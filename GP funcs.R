# Functions to sample from and predict values of a Gaussian process.

# Covariance function - squared exponential
K <- function(x1, x2, l, sigma_f) {
  sigma_f^2 * exp(-outer(x1, x2, function(x1, x2) (x1 - x2)^2) / (2 * l^2))
}

# Function to sample from a GP
samp_GP <- function(K, n_obs, mu = 0) {
  # Square-root of covariance matrix multiplied by standard normal vector.
  # Independent variance added for stability of Cholesky decomposition.
  # Transpose of Cholesky taken as upper triangular implementation rather than
  # lower.
  t(chol(K + 1e-10 * diag(n_obs))) %*% rnorm(n_obs) + mu
}

# Plot samples from a GP
plot_GP_samps = function(l, sigma_f, n_samps) {
  # Set number of samples
  n_obs = 100

  # Grid to sample/predict values over
  x <- seq(0, 1, len = n_obs)
  
  # Plot samples
  matplot(
    x, matrix(c(0, 1.96 * sigma_f, -1.96 * sigma_f), n_obs, 3, T),
    main = "Samples from a Gaussian process", 
    xlab = "x", ylab = "f(x)", t = 'l', col = 2, lty = c(1, 2, 2)
  )
  
  # Sample from GP and plot
  for (r in 1:n_samps) lines(x, samp_GP(K(x, x, l, sigma_f), n_obs))
}

# Function for Moore-Penrose pseudo-inverse via singular value decomposition
pinv <- function(mat) {
  svd_mat <- svd(mat)
  pinv_svs <- ifelse(svd_mat$d == 0, 0, 1 / svd_mat$d)
  svd_mat$v %*% diag(pinv_svs) %*% t(svd_mat$u)
}

# Function to predict from GP using pseudo inverse - deprecated
pred_GP_old <- function(x, x_star, l, sigma_f, sigma_n, n, y) {
  k_star <- K(x, x_star, l, sigma_f)
  k_inv <- pinv(K(x, x, l, sigma_f) + sigma_n^2 * diag(n))
  f_mean = t(k_star) %*% k_inv %*% y
  f_cov = K(x_star, x_star, l, sigma_f) - t(k_star) %*% k_inv %*% k_star
  list(f_mean, f_cov)
}

# Function to predict from GP using Cholesky decomposition
pred_GP <- function(x, x_star, l, sigma_f, sigma_n, n, y) {
  k_star <- K(x, x_star, l, sigma_f)
  L = t(chol(K(x, x, l, sigma_f) + sigma_n^2 * diag(n)))
  alpha = solve(t(L), solve(L, y))
  f_mean = t(k_star) %*% alpha
  v = solve(L, k_star)
  f_cov = K(x_star, x_star, l, sigma_f) - t(v) %*% v
  lml = t(y) %*% alpha / 2 + sum(diag(L)) - n/2 * log(2 * pi)
  list(f_mean, f_cov, lml)
}

# Trace of matrix-product
tr_mat_prod = function(mat1, mat2) {
  sum(sapply(1:n_obs, function(i) mat1[i, ] %*% mat2[, i]))
}

# Partial derivatives of log marginal likelihood w.r.t. hyperparameters
d_lml_d_hps = function(par, n, x, y) {
  l = par[1]
  sigma_f = par[2]
  sigma_n = par[3]
  K = K(x, x, l, sigma_f)
  K_y_inv = solve(K + sigma_n^2 * diag(n))
  alpha = K_y_inv %*% y
  mat_diff = alpha %*% t(alpha) - K_y_inv
  c(
    tr_mat_prod(mat_diff, sigma_f * K), # d_lml_d_sigma_f
    tr_mat_prod(mat_diff, K * outer(x, x, "-")^2 / l^3) / 2, # d_lml_d_l
    sum(diag(mat_diff)) / 2 # d_lml_d_sigma_n
  )
}

# Log marginal likelihood
lml <- function(par, n, x, y) {
  l = par[1]
  sigma_f = par[2]
  sigma_n = par[3]
  L = t(chol(K(x, x, l, sigma_f) + sigma_n^2 * diag(n)))
  alpha = solve(t(L), solve(L, y))
  t(y) %*% alpha / 2 + sum(diag(L)) - n/2 * log(2 * pi)
}

# Optimise log marginal likelihood w.r.t. hyperparameters
opt_lml = function() {
  opt = optim(par, lml, d_lml_d_hp, )
}

# Function to sample, predict, and plot
plot_GP_regression = function(n_obs, l, sigma_f, sigma_n) {
  # Number of data points to predict
  n_pred = 100
  
  # Grid to sample/predict values over
  x_obs <- seq(0, 1, len = n_obs)
  x_pred <- seq(0, 1, len = n_pred)
  
  # Sample from GP
  y_obs = samp_GP(K(x_obs, x_obs, l, sigma_f), n_obs) + 
    rnorm(n_obs, sd = sigma_n)
  
  # Predict mean and variance given samples
  y_pred = pred_GP(x_obs, x_pred, l, sigma_f, sigma_n, n_obs, y_obs)
  se = sqrt(diag(y_pred[[2]]))
  
  # Plot samples
  matplot(
    x_pred, matrix(y_pred[[1]], n_pred, 3) + cbind(0, 1.96 * se, -1.96 * se),
    main = "Gaussian process regression", xlab = "x", ylab = "f(x)", t = 'l',
    col = 2, lty = c(1, 2, 2)
  )
  points(x_obs, y_obs)
  
  # Sample from posterior GP and plot
  for (r in 1:3) lines(x_pred, samp_GP(y_pred[[2]], n_pred, y_pred[[1]]))
}
