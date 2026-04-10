# IV/GMM Estimation of Linear Parameters

Concentrates out linear parameters via instrumental variables.

## Usage

``` r
iv_estimate(X, Z, W, y, jacobian = NULL)
```

## Arguments

- X:

  Regressor matrix (N x K)

- Z:

  Instrument matrix (N x M)

- W:

  Weighting matrix (M x M)

- y:

  Dependent variable (N x 1)

- jacobian:

  Optional Jacobian of y w.r.t. nonlinear parameters (N x P)

## Value

List with parameters, residuals, covariances, and residual_jacobian
