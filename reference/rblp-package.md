# rblp: BLP Demand Estimation for Differentiated Products

Estimate, simulate, and analyze demand for differentiated products using
BLP-type random coefficients logit models (Berry, Levinsohn, and Pakes,
1995).

## Details

The rblp package provides tools for:

- Estimating random coefficients logit demand models via GMM

- Nested logit and random coefficients nested logit

- Joint demand and supply estimation

- Post-estimation analysis: elasticities, diversion ratios, markups

- Merger simulation and counterfactual analysis

- Consumer surplus computation

- Optimal instrument construction

- Micro moments (Conlon and Gortmaker, 2025)

- Data simulation for Monte Carlo studies

The main workflow is:

1.  Create formulations with
    [`blp_formulation`](https://hchulkim.github.io/rblp/reference/blp_formulation.md)

2.  Set up integration with
    [`blp_integration`](https://hchulkim.github.io/rblp/reference/blp_integration.md)

3.  Define the problem with
    [`blp_problem`](https://hchulkim.github.io/rblp/reference/blp_problem.md)

4.  Estimate with `problem$solve()`

5.  Analyze with result methods (e.g., `results$compute_elasticities()`)

Translated from the pyblp Python package by Conlon and Gortmaker (2020).

## References

Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in
Market Equilibrium. *Econometrica*, 63(4), 841-890.

Conlon, C., & Gortmaker, J. (2020). Best Practices for Differentiated
Products Demand Estimation with PyBLP. *RAND Journal of Economics*,
51(4), 1108-1161.

Nevo, A. (2000). A Practitioner's Guide to Estimation of
Random-Coefficients Logit Models of Demand. *Journal of Economics &
Management Strategy*, 9(4), 513-548.

## See also

Useful links:

- <https://github.com/hchulkim/rblp>

- Report bugs at <https://github.com/hchulkim/rblp/issues>

## Author

**Maintainer**: Hyoungchul Kim <hchulkim@virginia.edu>

Authors:

- Jeff Gortmaker (Author of original pyblp Python package)

- Chris Conlon (Co-author of pyblp methodology)
