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

      # Caching mechanism: R's optim() calls fn(theta) and gr(theta) as separate
      # function calls, even though they are evaluated at the same theta. In BLP,
      # each objective evaluation requires solving the full inner-loop contraction
      # mapping, which is expensive. This cache ensures that when optim requests
      # the gradient at the same theta it just used for the objective, we return
      # the already-computed gradient without re-solving the inner loop. The cache
      # uses an R environment (reference semantics) so it persists across calls.
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

      # Wrap the cached evaluator into separate fn/gr closures that optim expects.
      # Both closures share the same cache, so a call to fn followed by gr at the
      # same theta triggers only one underlying objective evaluation.
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

    # L-BFGS-B: limited-memory quasi-Newton with box constraints. This is the
    # default optimizer for BLP because it supports bounds on sigma parameters
    # (e.g., standard deviations must be non-negative) and uses gradient
    # information to approximate the Hessian via a limited-memory BFGS update.
    # factr = 0 disables the default convergence criterion based on relative
    # function change, relying instead on the projected gradient norm (pgtol).
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

    # BFGS: full-memory quasi-Newton without box constraints. More accurate
    # Hessian approximation than L-BFGS-B (stores the full inverse Hessian
    # rather than a limited number of update vectors), but cannot enforce
    # parameter bounds. Suitable when all nonlinear parameters are unrestricted.
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

    # Nelder-Mead (simplex): derivative-free optimizer. Does not use gradient
    # information, so it is robust to non-smooth or noisy objectives but
    # converges more slowly. Useful as a diagnostic or when analytic gradients
    # are unavailable. Higher default maxit because each simplex step is cheap
    # but many steps are needed for convergence.
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

    # nlminb: PORT trust-region optimizer with box constraints. An alternative
    # to L-BFGS-B that uses a trust-region approach rather than line search,
    # which can be more robust near constraint boundaries. Particularly useful
    # when parameters are near their bounds (e.g., sigma close to zero).
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
