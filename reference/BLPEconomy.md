# BLP Economy Base Class

Base class storing product/agent data and shared computations.

## Methods

### Public methods

- [`BLPEconomy$new()`](#method-BLPEconomy-new)

- [`BLPEconomy$get_market_data()`](#method-BLPEconomy-get_market_data)

- [`BLPEconomy$compute_logit_delta()`](#method-BLPEconomy-compute_logit_delta)

- [`BLPEconomy$print()`](#method-BLPEconomy-print)

- [`BLPEconomy$clone()`](#method-BLPEconomy-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    BLPEconomy$new(
      product_formulations,
      product_data,
      agent_formulation = NULL,
      agent_data = NULL,
      integration = NULL,
      rc_types = NULL,
      epsilon_scale = 1,
      costs_type = "linear",
      add_exogenous = TRUE
    )

------------------------------------------------------------------------

### Method `get_market_data()`

#### Usage

    BLPEconomy$get_market_data(market_id)

------------------------------------------------------------------------

### Method `compute_logit_delta()`

#### Usage

    BLPEconomy$compute_logit_delta(rho = NULL)

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

#### Usage

    BLPEconomy$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPEconomy$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
