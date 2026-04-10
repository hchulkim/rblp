# Mixtape Sessions Cereal Product Data

Simplified Nevo (2000) cereal data from the Mixtape Sessions demand
estimation exercises. Contains 2256 product-market observations with raw
quantities and a single price instrument.

## Usage

``` r
mixtape_products
```

## Format

A data frame with 2256 rows and the following columns:

- market:

  Market identifier (city-quarter)

- product:

  Product identifier (firm-brand)

- mushy:

  Mushiness indicator (0/1)

- servings_sold:

  Total servings sold in the market

- city_population:

  City population

- price_per_serving:

  Price per serving (dollars)

- price_instrument:

  Excluded instrument for price

## Source

<https://github.com/Mixtape-Sessions/Demand-Estimation>

## References

Conlon, C. & Gortmaker, J. (2020). Best Practices for Differentiated
Products Demand Estimation with pyblp. *RAND Journal of Economics*,
51(4), 1108-1161.
