# Functions to sample from and predict values of a Gaussian process.

# Covariance function - squared exponential
K <- function(x1, x2, l, sigma_f) {
  sigma_f^2 * exp(-outer(x1, x2, function(x1, x2) (x1 - x2)^2) / (2 * l^2))
}

# Function to sample from the GP
samp_GP <- function(x, l, sigma_f, n_samp) {
  # Square-root of covariance matrix multiplied by standard normal vector - add
  # independent variance for stability of Cholesky decomposition
  t(chol(K(x, x, l, sigma_f) + 1e-10 * diag(n_samp))) %*% rnorm(n_samp)
}

# Plot samples from a GP
plot_samp = function(n_samp, n_real, l, sigma_f) {
  # Grid to sample/predict values over
  x <- seq(0, 1, len = n_samp)
  
  # Plot samples
  matplot(
    x, matrix(c(0, 1.96 * sigma_f, -1.96 * sigma_f), n_samp, 3, T),
    main = "Samples from a Gaussian process", 
    xlab = "x", ylab = "f(x)", t = 'l', col = 2, lty = c(1, 2, 2)
  )
  
  # Sample from GP and plot
  for (r in 1:n_real) {
    y_grid = samp_GP(x, l, sigma_f, n_samp)
    lines(x, y_grid)
  }
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

# Function to sample, predict, and plot
plot_GP = function(n_samp, n_pred, l, sigma_f, sigma_n) {
  # Grid to sample/predict values over
  x_grid <- seq(0, 1, len = n_samp)
  x_pred <- seq(0, 1, len = n_pred)
  
  # Sample from GP
  y_grid = samp_GP(x_grid, l, sigma_f, n_samp)
  
  # Predict mean and variance given samples
  y_pred = pred_GP(x_grid, x_pred, l, sigma_f, sigma_n, n_samp, y_grid)
  
  # Plot samples
  plot(
    x_grid, y_grid, main = "Gaussian process regression", xlab = "x", 
    ylab = "f(x)"
  )
  lines(x_pred, y_pred[[1]], col = 2)
  lines(x_pred, y_pred[[1]] + 1.96 * sqrt(diag(y_pred[[2]])), col = 2, lty = 2)
  lines(x_pred, y_pred[[1]] - 1.96 * sqrt(diag(y_pred[[2]])), col = 2, lty = 2)
}
