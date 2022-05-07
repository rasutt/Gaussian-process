Gaussian Process
================
Robin Aldridge-Sutton
07/05/2022

``` r
# Functions to sample from and predict values of a Gaussian process.
source("GP funcs.R")
```

``` r
par(mfrow = c(2, 2))

for (i in 1:4)
  plot_GP(
    n_samp = 10, # Number of points to sample
    n_pred = 40, # Number of data points to predict
    l = 0.2, # Length scale
    sigma_f = 2, # Function standard deviation
    sigma_n = 0.1 # Noise standard deviation
  )
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
