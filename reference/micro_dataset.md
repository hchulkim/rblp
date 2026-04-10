# Create a Micro Dataset

Create a Micro Dataset

## Usage

``` r
micro_dataset(name, observations, compute_weights, market_ids = NULL)
```

## Arguments

- name:

  Dataset name

- observations:

  Number of observations

- compute_weights:

  Function(market_id, products, agents) -\> I x J weights

- market_ids:

  Optional market subset

## Value

A MicroDataset object
