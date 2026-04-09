#' @title Optimization Configuration
#' @description Configuration for the GMM outer-loop optimization.
#' @export
BLPOptimization <- R6::R6Class("BLPOptimization",
  public = list(
    #' @description Create optimization configuration
    #' @param method Character: optimizer name
    #' @param method_options Named list passed to the optimizer
    #' @param compute_gradient Logical: use analytic gradients?
    initialize = function(method = "l-bfgs-b", method_options = list(),
                          compute_gradient = TRUE) {
      valid <- c("l-bfgs-b", "bfgs", "nelder-mead", "nlminb", "return")
      if (!method %in% valid) {
        stop(sprintf("Invalid method '%s'. Choose from: %s",
                     method, paste(valid, collapse = ", ")))
      }
      private$method_ <- method
      private$options_ <- method_options
      private$compute_gradient_ <- compute_gradient
    },

    #' @description Run optimization
    #' @param initial Starting parameter values
    #' @param bounds List with lower and upper bound vectors
    #' @param objective_function Function: theta -> list(objective, gradient)
    #' @return List with values, converged, iterations, evaluations
    optimize = function(initial, bounds, objective_function) {
      if (private$method_ == "return") {
        return(list(values = initial, converged = TRUE,
                    iterations = 0L, evaluations = 0L))
      }

      # Cache to avoid double computation when optim calls fn and gr separately
      cache <- new.env(parent = emptyenv())
      cache$theta <- NULL
      cache$result <- NULL
      cache$evaluations <- 0L

      get_cached <- function(theta) {
        if (!identical(theta, cache$theta)) {
          cache$theta <- theta
          cache$result <- objective_function(theta)
          cache$evaluations <- cache$evaluations + 1L
          if (isTRUE(getOption("rblp.verbose", TRUE)) && cache$evaluations %% 25 == 0) {
            rblp_message(sprintf("  Optimization eval %d: objective = %.8e",
                                 cache$evaluations, cache$result$objective))
          }
        }
        cache$result
      }

      fn <- function(theta) get_cached(theta)$objective
      gr <- NULL
      if (private$compute_gradient_) {
        gr <- function(theta) get_cached(theta)$gradient
      }

      switch(private$method_,
        "l-bfgs-b" = private$optimize_lbfgsb(initial, bounds, fn, gr, cache),
        "bfgs"     = private$optimize_bfgs(initial, fn, gr, cache),
        "nelder-mead" = private$optimize_nm(initial, fn, cache),
        "nlminb"   = private$optimize_nlminb(initial, bounds, fn, gr, cache)
      )
    },

    #' @description Print optimization configuration
    print = function(...) {
      cat(sprintf("BLPOptimization: %s (gradient=%s)\n",
                  private$method_, private$compute_gradient_))
      invisible(self)
    }
  ),
  private = list(
    method_ = NULL,
    options_ = NULL,
    compute_gradient_ = TRUE,

    optimize_lbfgsb = function(initial, bounds, fn, gr, cache) {
      ctrl <- list(maxit = 1000L, factr = 0, pgtol = 1e-8)
      ctrl <- modifyList(ctrl, private$options_)

      result <- stats::optim(
        par = initial, fn = fn, gr = gr,
        method = "L-BFGS-B",
        lower = bounds$lower, upper = bounds$upper,
        control = ctrl
      )

      list(
        values = result$par,
        converged = result$convergence == 0,
        iterations = result$counts[1] %||% 0L,
        evaluations = cache$evaluations
      )
    },

    optimize_bfgs = function(initial, fn, gr, cache) {
      ctrl <- list(maxit = 1000L)
      ctrl <- modifyList(ctrl, private$options_)

      result <- stats::optim(
        par = initial, fn = fn, gr = gr,
        method = "BFGS", control = ctrl
      )

      list(
        values = result$par,
        converged = result$convergence == 0,
        iterations = result$counts[1] %||% 0L,
        evaluations = cache$evaluations
      )
    },

    optimize_nm = function(initial, fn, cache) {
      ctrl <- list(maxit = 5000L)
      ctrl <- modifyList(ctrl, private$options_)

      result <- stats::optim(
        par = initial, fn = fn,
        method = "Nelder-Mead", control = ctrl
      )

      list(
        values = result$par,
        converged = result$convergence == 0,
        iterations = result$counts[1] %||% 0L,
        evaluations = cache$evaluations
      )
    },

    optimize_nlminb = function(initial, bounds, fn, gr, cache) {
      ctrl <- list(iter.max = 1000L, eval.max = 2000L)
      ctrl <- modifyList(ctrl, private$options_)

      result <- stats::nlminb(
        start = initial, objective = fn, gradient = gr,
        lower = bounds$lower, upper = bounds$upper,
        control = ctrl
      )

      list(
        values = result$par,
        converged = result$convergence == 0,
        iterations = result$iterations %||% 0L,
        evaluations = cache$evaluations
      )
    }
  )
)

#' Create optimization configuration
#'
#' @param method Character: "l-bfgs-b" (default), "bfgs", "nelder-mead", "nlminb", or "return"
#' @param method_options Named list of options passed to the optimizer
#' @param compute_gradient Logical: use analytic gradients? (default TRUE)
#' @return A BLPOptimization object
#' @export
#' @examples
#' opt <- blp_optimization("l-bfgs-b")
blp_optimization <- function(method = "l-bfgs-b", method_options = list(),
                              compute_gradient = TRUE) {
  BLPOptimization$new(method, method_options, compute_gradient)
}
