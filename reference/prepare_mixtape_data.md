# Prepare Mixtape product data for rblp estimation

Transforms raw Mixtape product data into the format required by
[`blp_problem`](https://hchulkim.github.io/rblp/reference/blp_problem.md):
computes market shares from quantities, renames columns, and sets up
instruments.

## Usage

``` r
prepare_mixtape_data(products, servings_per_person = 90)
```

## Arguments

- products:

  Data frame from
  [`load_mixtape_products`](https://hchulkim.github.io/rblp/reference/load_mixtape_products.md)

- servings_per_person:

  Potential servings per person per quarter (default: 90, i.e. one
  serving per day for 90 days)

## Value

Data frame ready for
[`blp_problem`](https://hchulkim.github.io/rblp/reference/blp_problem.md)
