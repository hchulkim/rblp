# rblp Development Notes

## Package Overview
rblp is an R package implementing BLP (Berry, Levinsohn, Pakes 1995) demand estimation for differentiated products. It is translated from pyblp (Conlon & Gortmaker 2020).

## Key Architecture
- R6 classes: `BLPFormulation`, `BLPEconomy`, `BLPProblem`, `BLPResults`, `BLPSimulation`
- `BLPProblem` extends `BLPEconomy` (inherits data setup, adds solve methods)
- Main files: `R/problem.R` (solver), `R/economy.R` (base class), `R/results.R` (post-estimation), `R/formulation.R` (model specification), `R/simulation.R` (data generation)

## Testing
```bash
Rscript -e "testthat::test_dir('tests/testthat')"
Rscript inst/examples/nevo_example.R
```

## Key Conventions
- Product data must have: `market_ids`, `shares`, `prices`, `firm_ids`, `demand_instruments0..N`
- Agent data must have: `market_ids`, `weights`, plus demographic columns
- `blp_formulation(~ prices + x, absorb = ~ product_ids)` for fixed effects
- Integration: `blp_integration("product", size = N)` for Gauss-Hermite
- Optimization: `blp_optimization("l-bfgs-b", method_options = list(maxit = 200))`

## pyblp Correspondence
- pyblp `Formulation('0 + prices', absorb='C(product_ids)')` → `blp_formulation(~ 0 + prices, absorb = ~ product_ids)`
- pyblp `Problem([f1, f2, f3], product_data, agent_formulation, agent_data)` → `blp_problem(list(f1, f2), products, agent_formulation = f3, agent_data = agents)`
