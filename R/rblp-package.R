# Package-level documentation and imports. This file declares the roxygen
# preamble for the rblp package and defines shared utilities like %||%.

#' @title rblp: BLP Demand Estimation for Differentiated Products
#'
#' @description
#' Estimate, simulate, and analyze demand for differentiated products using
#' BLP-type random coefficients logit models (Berry, Levinsohn, and Pakes, 1995).
#'
#' @details
#' The rblp package provides tools for:
#' \itemize{
#'   \item Estimating random coefficients logit demand models via GMM
#'   \item Nested logit and random coefficients nested logit
#'   \item Joint demand and supply estimation
#'   \item Post-estimation analysis: elasticities, diversion ratios, markups
#'   \item Merger simulation and counterfactual analysis
#'   \item Consumer surplus computation
#'   \item Optimal instrument construction
#'   \item Micro moments (Conlon and Gortmaker, 2025)
#'   \item Data simulation for Monte Carlo studies
#' }
#'
#' The main workflow is:
#' \enumerate{
#'   \item Create formulations with \code{\link{blp_formulation}}
#'   \item Set up integration with \code{\link{blp_integration}}
#'   \item Define the problem with \code{\link{blp_problem}}
#'   \item Estimate with \code{problem$solve()}
#'   \item Analyze with result methods (e.g., \code{results$compute_elasticities()})
#' }
#'
#' Translated from the pyblp Python package by Conlon and Gortmaker (2020).
#'
#' @references
#' Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in Market
#' Equilibrium. \emph{Econometrica}, 63(4), 841-890.
#'
#' Conlon, C., & Gortmaker, J. (2020). Best Practices for Differentiated
#' Products Demand Estimation with PyBLP. \emph{RAND Journal of Economics},
#' 51(4), 1108-1161.
#'
#' Nevo, A. (2000). A Practitioner's Guide to Estimation of Random-Coefficients
#' Logit Models of Demand. \emph{Journal of Economics & Management Strategy},
#' 9(4), 513-548.
#'
#' @docType package
#' @name rblp-package
#' @aliases rblp
#' @importFrom R6 R6Class
#' @importFrom Matrix bdiag
#' @importFrom MASS ginv mvrnorm
#' @importFrom stats model.matrix optim nlminb qnorm rnorm runif pchisq
#'   setNames complete.cases sd var cov
#' @importFrom methods is
#' @importFrom parallel mclapply detectCores
#' @importFrom utils head tail
"_PACKAGE"

# Null-coalescing operator (internal utility, no roxygen to avoid Rd name issues)
`%||%` <- function(x, y) if (is.null(x)) y else x
