#' @title BLP Model Formulation
#' @description Configuration for designing model matrices from R formulas.
#' @export
BLPFormulation <- R6::R6Class("BLPFormulation",
  public = list(
    #' @description Create a new formulation
    #' @param formula Formula string or R formula
    #' @param absorb Optional formula for absorbed fixed effects
    initialize = function(formula, absorb = NULL) {
      if (is.character(formula)) {
        formula <- trimws(formula)
        if (!grepl("~", formula)) formula <- paste("~", formula)
        private$formula_ <- stats::as.formula(formula)
      } else {
        private$formula_ <- formula
      }
      if (!is.null(absorb)) {
        if (is.character(absorb)) {
          if (!grepl("~", absorb)) absorb <- paste("~", absorb)
          private$absorb_ <- stats::as.formula(absorb)
        } else {
          private$absorb_ <- absorb
        }
      }
    },

    #' @description Build design matrix from data
    #' @param data Data frame
    #' @return Matrix with named columns
    build_matrix = function(data) {
      # Convert the R formula into a numeric design matrix using R's standard
      # model.matrix machinery. This handles intercepts, factor expansion,
      # interactions, and transformations automatically from the formula syntax.
      mf <- stats::model.frame(private$formula_, data = data, na.action = stats::na.pass)
      mm <- stats::model.matrix(private$formula_, data = mf)

      # Track which columns involve prices or shares -- needed downstream to
      # identify endogenous variables that require instruments in the IV/GMM step.
      nms <- colnames(mm)
      private$names_ <- nms
      private$has_prices_ <- any(grepl("prices", nms, fixed = TRUE))
      private$has_shares_ <- any(grepl("shares", nms, fixed = TRUE))

      # Frisch-Waugh-Lovell (FWL) demeaning: absorbing fixed effects by
      # subtracting group means is algebraically equivalent to including a full
      # set of group dummies but avoids constructing the (potentially huge)
      # dummy matrix. For each absorbed factor variable, we compute the
      # within-group mean of every column and subtract it, projecting out the
      # group-level variation. After demeaning, the intercept is identically
      # zero and is dropped.
      if (!is.null(private$absorb_)) {
        fe_vars <- all.vars(private$absorb_)
        for (fv in fe_vars) {
          if (fv %in% names(data)) {
            grp <- data[[fv]]
            for (j in seq_len(ncol(mm))) {
              gm <- tapply(mm[, j], grp, mean)
              mm[, j] <- mm[, j] - gm[match(grp, names(gm))]
            }
          }
        }
        # The intercept is absorbed into the fixed effects, so remove it
        # to avoid a column of zeros that would cause rank deficiency.
        ic <- which(colnames(mm) == "(Intercept)")
        if (length(ic) > 0) mm <- mm[, -ic, drop = FALSE]
        private$names_ <- colnames(mm)
      }

      mm
    },

    #' @description Get column names
    #' @return Character vector
    get_names = function() private$names_,

    #' @description Check if prices are in the formulation
    has_prices = function() private$has_prices_,

    #' @description Check if shares are in the formulation
    has_shares = function() private$has_shares_,

    #' @description Get absorb formula (if any)
    get_absorb = function() private$absorb_,

    #' @description Print formulation
    print = function(...) {
      cat("BLPFormulation: ", deparse(private$formula_), "\n")
      if (!is.null(private$absorb_)) {
        cat("  Absorb: ", deparse(private$absorb_), "\n")
      }
      invisible(self)
    }
  ),
  private = list(
    formula_ = NULL,
    absorb_ = NULL,
    names_ = NULL,
    has_prices_ = FALSE,
    has_shares_ = FALSE
  )
)

#' Create a BLP formulation
#'
#' @param formula Formula string (e.g., "1 + prices + sugar") or R formula
#' @param absorb Optional formula for absorbed fixed effects (e.g., "~ product_id")
#' @return A BLPFormulation object
#' @export
#' @examples
#' f1 <- blp_formulation("1 + prices + sugar + mushy")
#' f2 <- blp_formulation("0 + prices", absorb = "~ product_ids")
blp_formulation <- function(formula, absorb = NULL) {
  BLPFormulation$new(formula, absorb)
}
