Gaussian Process
================
Robin Aldridge-Sutton
07/05/2022

``` r
# Functions to sample from and predict values of a Gaussian process.
source("GP funcs.R")

# GP parameters
true_l <- 0.2 # Length scale
true_sigma_f <- 2 # Function standard deviation
true_sigma_n <- 0.1 # Noise standard deviation
true_pars <- c(true_l, true_sigma_f, true_sigma_n)

# Number of data points to generate
n_grid <- 10
n_pred = 40

# Grid to sample/predict values over
x_grid <- seq(0, 1, len = n_grid)
x_pred <- seq(0, 1, len = n_pred)

# Sample from GP
y_grid = samp_GP(x_grid, true_l, true_sigma_f, n_grid, true_sigma_n)

# Predict mean and variance given samples
y_pred = pred_GP(x_grid, x_pred, true_l, true_sigma_f, true_sigma_n, n_grid, y_grid)

# Plot samples
plot(x_grid, y_grid)
lines(x_pred, y_pred[[1]], col = 2)
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
