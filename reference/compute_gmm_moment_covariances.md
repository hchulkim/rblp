# Compute GMM Moment Covariances

Compute GMM Moment Covariances

## Usage

``` r
compute_gmm_moment_covariances(
  u_list,
  Z_list,
  type = "robust",
  clustering_ids = NULL
)
```

## Arguments

- u_list:

  List of residual vectors (one per equation)

- Z_list:

  List of instrument matrices (one per equation)

- type:

  "robust", "clustered", or "unadjusted"

- clustering_ids:

  Optional clustering identifiers

## Value

Moment covariance matrix S
