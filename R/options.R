#' Get or Set rblp Package Options
#'
#' @param ... Named arguments to set options. If empty, returns all current options.
#' @return Invisibly returns the previous values of changed options, or all options if none changed.
#' @export
#' @examples
#' rblp_options(verbose = FALSE)
#' rblp_options()
rblp_options <- function(...) {
  args <- list(...)
  defaults <- list(
    rblp.digits = 7L,
    rblp.verbose = TRUE,
    rblp.verbose_output = "",
    rblp.pseudo_inverses = TRUE,
    rblp.collinear_atol = 1e-10,
    rblp.collinear_rtol = 1e-10,
    rblp.psd_atol = 1e-8,
    rblp.psd_rtol = 1e-8,
    rblp.finite_differences_epsilon = sqrt(.Machine$double.eps),
    rblp.weights_tol = 1e-10,
    rblp.micro_computation_chunks = 1L,
    rblp.num_processes = 1L
  )

  if (length(args) == 0) {
    current <- lapply(names(defaults), function(nm) {
      getOption(nm, defaults[[nm]])
    })
    names(current) <- sub("^rblp\\.", "", names(defaults))
    return(current)
  }

  old <- list()
  for (nm in names(args)) {
    full_nm <- paste0("rblp.", nm)
    old[[nm]] <- getOption(full_nm)
    opts <- list()
    opts[[full_nm]] <- args[[nm]]
    do.call(options, opts)
  }
  invisible(old)
}

.onLoad <- function(libname, pkgname) {
  opts <- list(
    rblp.digits = 7L,
    rblp.verbose = TRUE,
    rblp.pseudo_inverses = TRUE,
    rblp.collinear_atol = 1e-10,
    rblp.collinear_rtol = 1e-10,
    rblp.psd_atol = 1e-8,
    rblp.psd_rtol = 1e-8,
    rblp.finite_differences_epsilon = sqrt(.Machine$double.eps),
    rblp.weights_tol = 1e-10,
    rblp.micro_computation_chunks = 1L,
    rblp.num_processes = 1L
  )
  toset <- !(names(opts) %in% names(options()))
  if (any(toset)) options(opts[toset])
}
