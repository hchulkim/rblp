# Compute GMM Parameter Covariances

Compute GMM Parameter Covariances

## Usage

``` r
compute_gmm_parameter_covariances(W, S, G, se_type = "robust")
```

## Arguments

- W:

  Weighting matrix

- S:

  Moment covariance matrix

- G:

  Moment Jacobian (d_moments / d_theta)

- se_type:

  "robust", "clustered", or "unadjusted"

## Value

Parameter covariance matrix V
