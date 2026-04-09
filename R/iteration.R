#' @title Iteration Configuration
#' @description Configuration for the BLP contraction mapping fixed-point iteration.
#' @export
BLPIteration <- R6::R6Class("BLPIteration",
  public = list(
    #' @description Create iteration configuration
    #' @param method Character: "simple", "squarem", or "return"
    #' @param method_options Named list of options
    initialize = function(method = "squarem", method_options = list()) {
      valid <- c("simple", "squarem", "return")
      if (!method %in% valid) {
        stop(sprintf("Invalid method '%s'. Choose from: %s",
                     method, paste(valid, collapse = ", ")))
      }
      private$method_ <- method
      defaults <- list(
        max_evaluations = 5000L,
        atol = 1e-14,
        rtol = 0,
        scheme = 3L,
        step_min = 1,
        step_max = 1,
        step_factor = 4
      )
      private$options_ <- modifyList(defaults, method_options)
    },

    #' @description Run fixed-point iteration
    #' @param initial Starting values
    #' @param contraction Function: x -> list(x_new, weights) or just x_new
    #' @return List with values, converged, iterations, evaluations
    iterate = function(initial, contraction) {
      switch(private$method_,
        "simple"  = private$iterate_simple(initial, contraction),
        "squarem" = private$iterate_squarem(initial, contraction),
        "return"  = list(values = initial, converged = TRUE,
                         iterations = 0L, evaluations = 0L)
      )
    },

    #' @description Print iteration configuration
    print = function(...) {
      cat(sprintf("BLPIteration: %s (atol=%.1e, max_eval=%d)\n",
                  private$method_, private$options_$atol,
                  private$options_$max_evaluations))
      invisible(self)
    }
  ),
  private = list(
    method_ = NULL,
    options_ = NULL,

    iterate_simple = function(initial, contraction) {
      x <- initial
      opts <- private$options_
      converged <- FALSE

      for (k in seq_len(opts$max_evaluations)) {
        result <- contraction(x)
        x_new <- if (is.list(result)) result[[1]] else result
        diff <- max(abs(x_new - x))
        x_norm <- max(abs(x))
        tol <- opts$atol + opts$rtol * x_norm

        if (diff < tol) {
          converged <- TRUE
          x <- x_new
          break
        }
        x <- x_new
      }

      list(values = x, converged = converged,
           iterations = k, evaluations = k)
    },

    iterate_squarem = function(initial, contraction) {
      # SQUAREM acceleration (Varadhan & Roland 2008)
      x0 <- initial
      opts <- private$options_
      converged <- FALSE
      step_max <- opts$step_max
      evals <- 0L

      for (k in seq_len(opts$max_evaluations %/% 3L + 1L)) {
        # Step 1: x1 = T(x0)
        result1 <- contraction(x0)
        x1 <- if (is.list(result1)) result1[[1]] else result1
        evals <- evals + 1L

        r <- x1 - x0

        # Check convergence after first contraction
        diff <- max(abs(r))
        x_norm <- max(abs(x0))
        tol <- opts$atol + opts$rtol * x_norm
        if (diff < tol) {
          converged <- TRUE
          x0 <- x1
          break
        }

        # Step 2: x2 = T(x1)
        result2 <- contraction(x1)
        x2 <- if (is.list(result2)) result2[[1]] else result2
        evals <- evals + 1L

        v <- (x2 - x1) - r

        # Compute step length
        r_norm <- sqrt(sum(r^2))
        v_norm <- sqrt(sum(v^2))

        if (v_norm < .Machine$double.eps) {
          # No acceleration possible
          x0 <- x2
          next
        }

        alpha <- -r_norm / v_norm

        # Bound alpha
        alpha <- min(max(alpha, -step_max), -opts$step_min)

        # Extrapolate
        x_new <- x0 - 2 * alpha * r + alpha^2 * v

        # Step 3: stabilize with one more contraction
        result3 <- contraction(x_new)
        x_new <- if (is.list(result3)) result3[[1]] else result3
        evals <- evals + 1L

        # Check if we need to backstep
        diff_new <- max(abs(x_new - x0))
        diff_x2 <- max(abs(x2 - x0))

        if (!all(is.finite(x_new)) || diff_new > 2 * diff_x2) {
          # Backstep: use x2 instead
          x0 <- x2
        } else {
          x0 <- x_new
        }

        # Update step_max for scheme 3
        if (opts$scheme == 3 && abs(alpha) >= step_max) {
          step_max <- step_max * opts$step_factor
        }
      }

      list(values = x0, converged = converged,
           iterations = k, evaluations = evals)
    }
  )
)

#' Create iteration configuration
#'
#' @param method Character: "simple", "squarem", or "return"
#' @param method_options Named list. Key options:
#'   \describe{
#'     \item{max_evaluations}{Maximum contraction evaluations (default 5000)}
#'     \item{atol}{Absolute convergence tolerance (default 1e-14)}
#'     \item{rtol}{Relative convergence tolerance (default 0)}
#'   }
#' @return A BLPIteration object
#' @export
#' @examples
#' iter <- blp_iteration("squarem", list(atol = 1e-14))
blp_iteration <- function(method = "squarem", method_options = list()) {
  BLPIteration$new(method, method_options)
}
