#' Package setup (internal)
#'
#' @useDynLib hexsmoothR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  op <- options()
  defaults <- list(
    hexsmoothR.verbose = TRUE,
    hexsmoothR.chunk_size = 1000L,
    hexsmoothR.use_sparse = TRUE
  )
  toset <- !(names(defaults) %in% names(op))
  if (any(toset)) options(defaults[toset])

  required_packages <- c("sf", "terra", "exactextractr")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    warning(
      "hexsmoothR: Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible()
}

#' Package startup message (internal)
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    packageStartupMessage(
      "hexsmoothR: 'Matrix' not available; sparse-matrix optimisation disabled."
    )
  }
}

#' Internal verbose-print helper
#'
#' Honours `options(hexsmoothR.verbose = ...)`. Use this instead of bare `cat()`
#' so users can silence the package globally.
#'
#' @param ... Passed to `cat()`.
#' @keywords internal
hex_msg <- function(...) {
  if (isTRUE(getOption("hexsmoothR.verbose", TRUE))) {
    cat(..., sep = "")
  }
  invisible()
}
