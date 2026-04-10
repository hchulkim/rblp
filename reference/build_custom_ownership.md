# Build Custom Ownership Matrix

Creates ownership matrices with custom specifications.

## Usage

``` r
build_custom_ownership(product_data, kappa = "standard")
```

## Arguments

- product_data:

  Data frame with market_ids and firm_ids

- kappa:

  "standard" (same firm = 1), "monopoly" (all 1), or "single" (identity)

## Value

Block-diagonal ownership matrix (N x N)
