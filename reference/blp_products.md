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

  Excluded demand-side instrument 0

- demand_instruments1:

  Excluded demand-side instrument 1

- demand_instruments2:

  Excluded demand-side instrument 2

- demand_instruments3:

  Excluded demand-side instrument 3

- demand_instruments4:

  Excluded demand-side instrument 4

- demand_instruments5:

  Excluded demand-side instrument 5

- demand_instruments6:

  Excluded demand-side instrument 6

- demand_instruments7:

  Excluded demand-side instrument 7

- supply_instruments0:

  Excluded supply-side instrument 0

- supply_instruments1:

  Excluded supply-side instrument 1

- supply_instruments2:

  Excluded supply-side instrument 2

- supply_instruments3:

  Excluded supply-side instrument 3

- supply_instruments4:

  Excluded supply-side instrument 4

- supply_instruments5:

  Excluded supply-side instrument 5

- supply_instruments6:

  Excluded supply-side instrument 6

- supply_instruments7:

  Excluded supply-side instrument 7

- supply_instruments8:

  Excluded supply-side instrument 8

- supply_instruments9:

  Excluded supply-side instrument 9

- supply_instruments10:

  Excluded supply-side instrument 10

- supply_instruments11:

  Excluded supply-side instrument 11

## Source

<https://pyblp.readthedocs.io/en/stable/>

## References

Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in
Market Equilibrium. *Econometrica*, 63(4), 841-890.
