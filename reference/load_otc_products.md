# OTC Drug Product Data (baby_BLP)

Product-level data from Lei Ma's baby_BLP project. Contains weekly sales
of 11 over-the-counter pain relief products across 2 stores and 48
weeks, yielding 96 store-week markets and 1,056 observations.

## Usage

``` r
load_otc_products()
```

## Format

A data frame with 1056 rows and the following columns:

- market_ids:

  Market identifier (store-week)

- firm_ids:

  Product identifier (1–11)

- shares:

  Market share (sales / market size)

- prices:

  Retail price

- promotion:

  Promotion indicator (0/1)

- cost:

  Wholesale cost

- product:

  Product factor (for fixed effects)

## Source

<https://github.com/leima0521/baby_BLP>

## References

Lei Ma, baby_BLP: A pedagogical implementation of BLP demand estimation.

## See also

[`vignette("baby-blp-replication")`](https://hchulkim.github.io/rblp/articles/baby-blp-replication.md)
for a walkthrough
