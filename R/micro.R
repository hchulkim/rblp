#' @title Micro Dataset
#' @description Defines a micro-level dataset for micro moments.
#' @export
#
# Implements the micro moments framework of Conlon & Gortmaker (2020).
# Traditional BLP uses only aggregate market shares as data; micro moments
# add individual-level survey data (e.g., Nielsen panel, consumer surveys)
# as additional moment conditions. This sharpens identification of the
# random coefficient distribution by leveraging information about WHO buys
# WHAT -- not just how much of each product is sold.
#
# The framework has three layers:
#   MicroDataset: defines the survey population (who was sampled, how many)
#   MicroPart: computes a weighted expectation (e.g., E[income | chose j])
#   MicroMoment: matches observed survey statistic to model prediction
MicroDataset <- R6::R6Class("MicroDataset",
  public = list(
    #' @field name Dataset name
    name = NULL,
    #' @field observations Number of observations in micro data
    observations = NULL,
    #' @field compute_weights Function(market_id, products, agents) -> I x J weights
    compute_weights = NULL,
    #' @field market_ids Optional subset of markets
    market_ids = NULL,

    #' @description Create micro dataset
    #' @param name Dataset name
    #' @param observations Number of micro observations
    #' @param compute_weights Function returning agent x product weight matrix
    #' @param market_ids Optional market subset
    #
    # A MicroDataset defines the sampling frame of a micro-level survey.
    # compute_weights returns an I x J matrix where entry (i,j) is the
    # probability that individual i purchasing product j would appear in
    # this survey. For a simple random sample of all buyers, weights are
    # all 1. For a survey restricted to certain products or demographics,
    # the weights encode the selection. The number of observations controls
    # the variance of the micro moment in the GMM weighting matrix.
    initialize = function(name, observations, compute_weights, market_ids = NULL) {
      self$name <- name
      self$observations <- observations
      self$compute_weights <- compute_weights
      self$market_ids <- market_ids
    },

    #' @description Print the dataset
    print = function(...) {
      cat(sprintf("MicroDataset: '%s' (%d observations)\n",
                  self$name, self$observations))
      invisible(self)
    }
  )
)

#' @title Micro Part
#' @description A single component of a micro moment calculation.
#' @export
#
# A MicroPart computes the model's prediction of a conditional expectation
# from the micro data. For example, "what is the average income of consumers
# who purchase cereal brand j?" The model predicts this via:
#   E[d_i | chose j] = sum_i w_i * P_ij * d_i / sum_i w_i * P_ij
# where P_ij is the choice probability from the BLP model, w_i are agent
# weights, and d_i is the demographic value. The compute_values function
# returns the I x J matrix of "values" (e.g., incomes), which gets
# weighted by choice probabilities and survey weights to form the
# model-predicted conditional moment.
MicroPart <- R6::R6Class("MicroPart",
  public = list(
    #' @field name Part name
    name = NULL,
    #' @field dataset MicroDataset object
    dataset = NULL,
    #' @field compute_values Function(market_id, products, agents) -> I x J values
    compute_values = NULL,

    #' @description Create a micro part
    #' @param name Part name
    #' @param dataset MicroDataset object
    #' @param compute_values Function returning agent x product value matrix
    initialize = function(name, dataset, compute_values) {
      self$name <- name
      self$dataset <- dataset
      self$compute_values <- compute_values
    },

    #' @description Compute aggregated part value across markets
    #' @param economy BLPEconomy object
    #' @param delta Mean utilities
    #' @param sigma Sigma matrix
    #' @param pi Pi matrix
    #' @param rho Rho vector
    #' @return Scalar aggregated value
    #
    # Evaluates the model-predicted conditional expectation across all
    # relevant markets. For each market, computes individual choice
    # probabilities P_ij from the current parameters, then forms the
    # weighted average:
    #   part_value = sum_{t,i,j} w_i * survey_wt(i,j) * P_ij * value(i,j)
    #              / sum_{t,i,j} w_i * survey_wt(i,j) * P_ij
    # This is the model's prediction of the micro statistic (e.g.,
    # "average income of brand j buyers") given the current parameter
    # values. The denominator normalizes by total predicted survey mass.
    compute = function(economy, delta, sigma = NULL, pi = NULL, rho = NULL) {
      market_list <- self$dataset$market_ids %||% economy$unique_market_ids
      numerator <- 0
      denominator <- 0

      for (t in market_list) {
        md <- economy$get_market_data(t)
        idx <- md$indices

        mkt <- BLPMarket$new(
          products = md$products,
          agents = md$agents,
          sigma = sigma, pi = pi, rho = rho,
          rc_types = economy$rc_types,
          epsilon_scale = economy$epsilon_scale
        )

        # Compute individual choice probabilities P (J x I) from the current
        # theta, delta, and mu. These probabilities are the core link between
        # the structural model and the micro data.
        mu <- mkt$compute_mu()
        prob_result <- mkt$compute_probabilities(delta[idx], mu)
        P <- if (is.list(prob_result)) prob_result$probabilities else prob_result

        # Survey weights (which consumer-product pairs appear in the data)
        # and value function (the statistic being averaged, e.g., income).
        weights_matrix <- self$dataset$compute_weights(t, md$products, md$agents)
        values_matrix <- self$compute_values(t, md$products, md$agents)

        # Accumulate the weighted conditional expectation. Each (i,j) pair
        # contributes: agent_weight * survey_weight * choice_prob * value.
        # The ratio numerator/denominator gives E[value | in survey, chose].
        w <- md$agents$weights
        for (i in seq_len(mkt$I)) {
          for (j in seq_len(mkt$J)) {
            wt <- weights_matrix[i, j] * P[j, i] * w[i]
            numerator <- numerator + wt * values_matrix[i, j]
            denominator <- denominator + wt
          }
        }
      }

      if (abs(denominator) > 1e-300) numerator / denominator else 0
    },

    #' @description Print the part
    print = function(...) {
      cat(sprintf("MicroPart: '%s'\n", self$name))
      invisible(self)
    }
  )
)

#' @title Micro Moment
#' @description Matches a micro data target to a model prediction.
#' @export
#
# A MicroMoment defines a single moment condition of the form:
#   g_micro(theta) = f(simulated_parts) - observed_value
# where observed_value is computed from survey micro data (e.g., "average
# income of Coke buyers = $45,000") and simulated_parts is the model's
# prediction of that same statistic given parameter theta. In GMM, this
# moment is stacked alongside the standard aggregate BLP moments
# E[Z' xi] = 0 to form a joint system. Adding micro moments tightens
# identification because they discipline the distribution of heterogeneity
# (who buys what) rather than just aggregate substitution patterns.
MicroMoment <- R6::R6Class("MicroMoment",
  public = list(
    #' @field name Moment name
    name = NULL,
    #' @field value Observed (target) value from micro data
    value = NULL,
    #' @field parts List of MicroPart objects
    parts = NULL,
    #' @field compute_value_fn Function of part values -> scalar
    compute_value_fn = NULL,
    #' @field compute_gradient_fn Function of part values -> gradient vector
    compute_gradient_fn = NULL,

    #' @description Create micro moment
    #' @param name Moment name
    #' @param value Observed scalar target
    #' @param parts Single MicroPart or list of MicroPart objects
    #' @param compute_value Function of part values -> scalar (default: identity)
    #' @param compute_gradient Function of part values -> gradient (default: 1)
    #
    # The moment condition is: g = compute_value(part_values) - value.
    # For simple moments (e.g., conditional mean), compute_value is identity
    # and there is one part. For ratio moments (e.g., choice probability
    # conditional on demographics), compute_value can be a ratio of two
    # parts. The compute_gradient function provides the Jacobian of
    # compute_value with respect to part values, needed for GMM standard
    # errors via the delta method.
    initialize = function(name, value, parts,
                          compute_value = NULL,
                          compute_gradient = NULL) {
      self$name <- name
      self$value <- value
      if (inherits(parts, "MicroPart")) parts <- list(parts)
      self$parts <- parts
      self$compute_value_fn <- compute_value %||% function(v) {
        if (length(v) == 1) v[[1]] else v[[1]]
      }
      self$compute_gradient_fn <- compute_gradient %||% function(v) rep(1, length(v))
    },

    #' @description Compute simulated moment value
    #' @param economy BLPEconomy object
    #' @param delta Mean utilities
    #' @param sigma Sigma matrix
    #' @param pi Pi matrix
    #' @param rho Rho vector
    #' @return Scalar simulated value
    #
    # Evaluates the model's prediction of this micro statistic at the
    # current parameter values. Each MicroPart is computed (producing a
    # conditional expectation), then compute_value_fn combines them into
    # a scalar. The GMM residual for this moment is:
    #   g = compute_simulated_value(...) - self$value
    # At the true parameters, g should be approximately zero.
    compute_simulated_value = function(economy, delta,
                                        sigma = NULL, pi = NULL, rho = NULL) {
      part_values <- lapply(self$parts, function(p) {
        p$compute(economy, delta, sigma, pi, rho)
      })
      self$compute_value_fn(part_values)
    },

    #' @description Print the moment
    print = function(...) {
      cat(sprintf("MicroMoment: '%s' (target = %.6f, %d parts)\n",
                  self$name, self$value, length(self$parts)))
      invisible(self)
    }
  )
)

#' Create a Micro Dataset
#'
#' @param name Dataset name
#' @param observations Number of observations
#' @param compute_weights Function(market_id, products, agents) -> I x J weights
#' @param market_ids Optional market subset
#' @return A MicroDataset object
#' @export
micro_dataset <- function(name, observations, compute_weights, market_ids = NULL) {
  MicroDataset$new(name, observations, compute_weights, market_ids)
}

#' Create a Micro Part
#'
#' @param name Part name
#' @param dataset MicroDataset object
#' @param compute_values Function(market_id, products, agents) -> I x J values
#' @return A MicroPart object
#' @export
micro_part <- function(name, dataset, compute_values) {
  MicroPart$new(name, dataset, compute_values)
}

#' Create a Micro Moment
#'
#' @param name Moment name
#' @param value Observed target value
#' @param parts Single MicroPart or list of MicroPart objects
#' @param compute_value Function of part values -> scalar
#' @param compute_gradient Function of part values -> gradient
#' @return A MicroMoment object
#' @export
micro_moment <- function(name, value, parts,
                          compute_value = NULL, compute_gradient = NULL) {
  MicroMoment$new(name, value, parts, compute_value, compute_gradient)
}
