# Nevo (2000) Cereal Agent Data

Agent-level data (simulated consumers) for the Nevo (2000) cereal study.
Contains integration nodes, weights, and demographic variables for 20
agents per market.

## Usage

``` r
nevo_agents
```

## Format

A data frame with 1880 rows and the following columns:

- market_ids:

  Market identifier (must match nevo_products)

- city_ids:

  City identifier

- quarter:

  Quarter

- weights:

  Integration weights (sum to 1 within each market)

- nodes0:

  Integration node dimension 0

- nodes1:

  Integration node dimension 1

- nodes2:

  Integration node dimension 2

- nodes3:

  Integration node dimension 3

- income:

  Log income

- income_squared:

  Squared log income

- age:

  Age (transformed)

- child:

  Child presence indicator (transformed)

## Source

<https://pyblp.readthedocs.io/en/stable/>

## References

Nevo, A. (2000). A Practitioner's Guide to Estimation of
Random-Coefficients Logit Models of Demand. *Journal of Economics &
Management Strategy*, 9(4), 513-548.
