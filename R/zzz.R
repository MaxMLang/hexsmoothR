#' Package setup (internal)
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  options(hexsmoothR.verbose = TRUE)
  options(hexsmoothR.chunk_size = 1000)
  options(hexsmoothR.use_sparse = TRUE)
  required_packages <- c("sf", "terra", "exactextractr", "data.table")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) warning("hexsmoothR: Missing required packages: ", paste(missing_packages, collapse = ", "))
  
  # Load the C++ library if available
  tryCatch({
    library.dynam("hexsmoothR", pkgname, libname)
    options(hexsmoothR.cpp_available = TRUE)
  }, error = function(e) {
    options(hexsmoothR.cpp_available = FALSE)
  })
}

#' Package startup message (internal)
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    packageStartupMessage("hexsmoothR: Matrix package not available. Sparse matrix optimization will be disabled.")
  }
  
  # Check C++ availability and show appropriate message
  if (!getOption("hexsmoothR.cpp_available", FALSE)) {
    packageStartupMessage("hexsmoothR: C++ library not available, will use R fallback implementations")
  }
} 