# BLP Market Computations

Market-level computations for BLP demand estimation. All probability,
share, delta, markup, and Jacobian computations happen here.

## Methods

### Public methods

- [`BLPMarket$new()`](#method-BLPMarket-new)

- [`BLPMarket$compute_random_coefficients()`](#method-BLPMarket-compute_random_coefficients)

- [`BLPMarket$compute_mu()`](#method-BLPMarket-compute_mu)

- [`BLPMarket$compute_probabilities()`](#method-BLPMarket-compute_probabilities)

- [`BLPMarket$compute_shares()`](#method-BLPMarket-compute_shares)

- [`BLPMarket$compute_delta()`](#method-BLPMarket-compute_delta)

- [`BLPMarket$compute_shares_by_xi_jacobian()`](#method-BLPMarket-compute_shares_by_xi_jacobian)

- [`BLPMarket$compute_shares_by_theta_jacobian()`](#method-BLPMarket-compute_shares_by_theta_jacobian)

- [`BLPMarket$compute_xi_by_theta_jacobian()`](#method-BLPMarket-compute_xi_by_theta_jacobian)

- [`BLPMarket$compute_eta()`](#method-BLPMarket-compute_eta)

- [`BLPMarket$compute_costs()`](#method-BLPMarket-compute_costs)

- [`BLPMarket$compute_equilibrium_prices()`](#method-BLPMarket-compute_equilibrium_prices)

- [`BLPMarket$compute_elasticities()`](#method-BLPMarket-compute_elasticities)

- [`BLPMarket$compute_diversion_ratios()`](#method-BLPMarket-compute_diversion_ratios)

- [`BLPMarket$compute_consumer_surplus()`](#method-BLPMarket-compute_consumer_surplus)

- [`BLPMarket$compute_hhi()`](#method-BLPMarket-compute_hhi)

- [`BLPMarket$compute_markups()`](#method-BLPMarket-compute_markups)

- [`BLPMarket$compute_profits()`](#method-BLPMarket-compute_profits)

- [`BLPMarket$clone()`](#method-BLPMarket-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    BLPMarket$new(
      products,
      agents,
      sigma = NULL,
      pi = NULL,
      rho = NULL,
      rc_types = NULL,
      epsilon_scale = 1,
      costs_type = "linear"
    )

------------------------------------------------------------------------

### Method `compute_random_coefficients()`

#### Usage

    BLPMarket$compute_random_coefficients()

------------------------------------------------------------------------

### Method `compute_mu()`

#### Usage

    BLPMarket$compute_mu(coefficients = NULL)

------------------------------------------------------------------------

### Method `compute_probabilities()`

#### Usage

    BLPMarket$compute_probabilities(delta, mu = NULL)

------------------------------------------------------------------------

### Method `compute_shares()`

#### Usage

    BLPMarket$compute_shares(probabilities, weights = NULL)

------------------------------------------------------------------------

### Method `compute_delta()`

#### Usage

    BLPMarket$compute_delta(initial_delta, iteration, fp_type = "safe_linear")

------------------------------------------------------------------------

### Method `compute_shares_by_xi_jacobian()`

#### Usage

    BLPMarket$compute_shares_by_xi_jacobian(probabilities, conditionals = NULL)

------------------------------------------------------------------------

### Method `compute_shares_by_theta_jacobian()`

#### Usage

    BLPMarket$compute_shares_by_theta_jacobian(
      probabilities,
      conditionals = NULL,
      sigma,
      pi_mat = NULL,
      rho = NULL
    )

------------------------------------------------------------------------

### Method `compute_xi_by_theta_jacobian()`

#### Usage

    BLPMarket$compute_xi_by_theta_jacobian(probabilities, conditionals = NULL)

------------------------------------------------------------------------

### Method `compute_eta()`

#### Usage

    BLPMarket$compute_eta(probabilities, ownership = NULL)

------------------------------------------------------------------------

### Method `compute_costs()`

#### Usage

    BLPMarket$compute_costs(probabilities, prices, ownership = NULL)

------------------------------------------------------------------------

### Method `compute_equilibrium_prices()`

#### Usage

    BLPMarket$compute_equilibrium_prices(
      costs,
      iteration,
      ownership = NULL,
      prices_init = NULL
    )

------------------------------------------------------------------------

### Method `compute_elasticities()`

#### Usage

    BLPMarket$compute_elasticities(probabilities, prices = NULL, name = "prices")

------------------------------------------------------------------------

### Method `compute_diversion_ratios()`

#### Usage

    BLPMarket$compute_diversion_ratios(probabilities, name = "prices")

------------------------------------------------------------------------

### Method `compute_consumer_surplus()`

#### Usage

    BLPMarket$compute_consumer_surplus(probabilities, delta, mu)

------------------------------------------------------------------------

### Method `compute_hhi()`

#### Usage

    BLPMarket$compute_hhi(shares = NULL, firm_ids = NULL)

------------------------------------------------------------------------

### Method `compute_markups()`

#### Usage

    BLPMarket$compute_markups(prices, costs)

------------------------------------------------------------------------

### Method `compute_profits()`

#### Usage

    BLPMarket$compute_profits(prices, shares, costs)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPMarket$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
