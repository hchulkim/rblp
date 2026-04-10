# BLP (1995) Automobile Product Data

Product-level data from the Berry, Levinsohn, and Pakes (1995) study of
the US automobile market. Contains 2217 product-market observations
across 20 annual markets (1971-1990).

## Usage

``` r
blp_products
```

## Format

A data frame with 2217 rows and the following columns:

- market_ids:

  Market year

- clustering_ids:

  Clustering identifier for standard errors

- car_ids:

  Vehicle identifier

- firm_ids:

  Manufacturer identifier

- region:

  Region of origin (US, Japan, Europe)

- shares:

  Market share

- prices:

  Price (10,000s of 1983 dollars)

- hpwt:

  Horsepower/weight ratio

- air:

  Air conditioning indicator

- mpd:

  Miles per dollar

- mpg:

  Miles per gallon

- space:

  Interior space (length x width)

- trend:

  Time trend

- demand_instruments0:

  Excluded demand-side instruments (8 total: 0-7)

- supply_instruments0:

  Excluded supply-side instruments (12 total: 0-11)

## Source

<https://pyblp.readthedocs.io/en/stable/>

## References

Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in
Market Equilibrium. *Econometrica*, 63(4), 841-890.
