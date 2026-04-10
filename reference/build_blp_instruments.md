# Build BLP Instruments (Sums of Characteristics)

Constructs excluded instruments using the traditional BLP approach: sums
of characteristics of other own-firm products and rival products.

## Usage

``` r
build_blp_instruments(X, market_ids, firm_ids)
```

## Arguments

- X:

  Matrix of exogenous product characteristics (N x K)

- market_ids:

  Market identifiers (length N)

- firm_ids:

  Firm identifiers (length N)

## Value

Matrix of instruments (N x 2K): own-firm and rival columns

## Examples

``` r
if (FALSE) { # \dontrun{
iv <- build_blp_instruments(X_exog, product_data$market_ids, product_data$firm_ids)
} # }
```
