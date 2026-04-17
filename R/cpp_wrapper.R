#' R fallback for N-order spatial smoothing
#'
#' Pure-R implementation matching the C++ algorithm in
#' `process_district_all_vars_n_orders()`. Used automatically by
#' [`smooth_variables()`] when the compiled library is unavailable.
#'
#' @param variable_values List of numeric vectors (one per variable).
#' @param neighbors List of length `n_orders`; each element is a list of
#'   length `n_cells` of integer neighbour indices.
#' @param weights List with `center_weight` and `neighbor_weights`.
#' @param hex_indices Integer vector of cell indices (1-based) to process.
#' @param var_names Character vector of variable names.
#' @return Named list with one entry per variable, each containing
#'   `raw` (the centre cell value), `weighted_combined` (the smoothed value),
#'   and `neighbors_<N><suffix>` (per-order neighbour means, e.g.
#'   `neighbors_1st`, `neighbors_2nd`, ...).
#' @keywords internal
#' @export
smooth_variables_r_fallback <- function(variable_values, neighbors, weights,
                                        hex_indices, var_names) {
  validate_smoothing_inputs(variable_values, neighbors, weights, var_names)

  n_vars     <- length(variable_values)
  n_process  <- length(hex_indices)
  n_orders   <- length(neighbors)
  center_w   <- weights$center_weight
  nb_weights <- weights$neighbor_weights

  result <- vector("list", n_vars)
  names(result) <- var_names

  for (v in seq_len(n_vars)) {
    values <- variable_values[[v]]

    raw_vals          <- rep(NA_real_, n_process)
    weighted_combined <- rep(NA_real_, n_process)
    order_means       <- vector("list", n_orders)
    for (order in seq_len(n_orders)) {
      order_means[[order]] <- rep(NA_real_, n_process)
    }

    for (i in seq_len(n_process)) {
      hex_idx <- hex_indices[i]
      if (hex_idx < 1 || hex_idx > length(neighbors[[1]])) next

      raw_vals[i] <- values[hex_idx]

      for (order in seq_len(n_orders)) {
        neighs <- neighbors[[order]][[hex_idx]]
        if (length(neighs) > 0) {
          order_vals <- values[neighs]
          if (any(!is.na(order_vals))) {
            order_means[[order]][i] <- mean(order_vals, na.rm = TRUE)
          }
        }
      }

      weighted_sum <- 0
      weight_sum   <- 0

      if (!is.na(values[hex_idx])) {
        weighted_sum <- weighted_sum + values[hex_idx] * center_w
        weight_sum   <- weight_sum + center_w
      }

      for (order in seq_len(n_orders)) {
        neighs <- neighbors[[order]][[hex_idx]]
        if (length(neighs) > 0) {
          ok <- !is.na(values[neighs])
          if (any(ok)) {
            n_ok <- sum(ok)
            weighted_sum <- weighted_sum + sum(values[neighs][ok]) * nb_weights[order]
            weight_sum   <- weight_sum + n_ok * nb_weights[order]
          }
        }
      }

      if (weight_sum > 0) {
        weighted_combined[i] <- weighted_sum / weight_sum
      }
    }

    var_result <- list(
      raw = raw_vals,
      weighted_combined = weighted_combined
    )
    for (order in seq_len(n_orders)) {
      suffix <- ordinal_suffix(order)
      var_result[[paste0("neighbors_", order, suffix)]] <- order_means[[order]]
    }

    result[[v]] <- var_result
  }

  result
}

#' Apply spatial smoothing to variables on a hexagonal grid
#'
#' Smooths variables using the topology produced by [`compute_topology()`].
#' Uses the compiled C++ implementation when available; if the C++ call
#' fails for an unexpected reason, falls back to the pure-R implementation
#' and re-raises clear validation errors.
#'
#' @param variable_values Named list of numeric vectors (one per variable).
#' @param neighbors List of neighbour lists (per order). Either the `neighbors`
#'   element of a topology, or the topology object itself.
#' @param weights List with `center_weight` and `neighbor_weights`. Either the
#'   `weights` element of a topology, or the topology object itself.
#' @param hex_indices Integer vector of cell indices to process. Defaults to
#'   all cells.
#' @param var_names Character vector of variable names. Defaults to
#'   `names(variable_values)`.
#' @return Named list with one entry per variable. Each entry contains:
#'   \itemize{
#'     \item `raw`: the original (unsmoothed) centre-cell values
#'     \item `weighted_combined`: weighted average of centre + all neighbours
#'     \item `neighbors_<N><suffix>`: mean of neighbours at order N
#'       (e.g. `neighbors_1st`, `neighbors_2nd`, `neighbors_3rd`)
#'   }
#' @export
smooth_variables <- function(variable_values, neighbors, weights,
                             hex_indices = NULL, var_names = NULL) {

  if (is.null(neighbors) || is.null(weights)) {
    stop("Both 'neighbors' and 'weights' must be provided", call. = FALSE)
  }

  if (is.list(weights) && "weights" %in% names(weights) &&
      !all(c("center_weight", "neighbor_weights") %in% names(weights))) {
    weights <- weights$weights
    hex_msg("Note: extracted weights from topology object\n")
  }
  if (is.list(neighbors) && "neighbors" %in% names(neighbors) &&
      !all(vapply(neighbors, is.list, logical(1)))) {
    neighbors <- neighbors$neighbors
    hex_msg("Note: extracted neighbors from topology object\n")
  }

  if (is.null(var_names)) var_names <- names(variable_values)
  if (is.null(var_names) || any(!nzchar(var_names))) {
    stop("var_names must be supplied (or variable_values must be named)", call. = FALSE)
  }

  validate_smoothing_inputs(variable_values, neighbors, weights, var_names)

  if (is.null(hex_indices)) {
    hex_indices <- seq_along(neighbors[[1]])
  }
  hex_indices <- as.integer(hex_indices)

  cpp_available <- isTRUE(getOption("hexsmoothR.cpp_available", TRUE)) &&
    exists("process_district_all_vars_n_orders",
           envir = asNamespace("hexsmoothR"),
           inherits = FALSE)

  if (cpp_available) {
    out <- tryCatch(
      process_district_all_vars_n_orders(
        variable_values, neighbors, weights, hex_indices, var_names
      ),
      error = function(e) {
        warning(
          "C++ smoothing failed (", e$message, "); using R fallback.",
          call. = FALSE
        )
        NULL
      }
    )
    if (!is.null(out)) return(out)
  }

  smooth_variables_r_fallback(
    variable_values, neighbors, weights, hex_indices, var_names
  )
}

#' Validate inputs to the smoothing functions (internal)
#' @keywords internal
validate_smoothing_inputs <- function(variable_values, neighbors, weights, var_names) {
  if (!is.list(variable_values) || length(variable_values) == 0) {
    stop("variable_values must be a non-empty list", call. = FALSE)
  }
  if (!all(vapply(variable_values, is.numeric, logical(1)))) {
    stop("variable_values must contain numeric vectors only", call. = FALSE)
  }
  sizes <- vapply(variable_values, length, integer(1))
  if (length(unique(sizes)) > 1) {
    stop(
      "All variables must be the same length. Got: ",
      paste(sizes, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.list(weights) ||
      !all(c("center_weight", "neighbor_weights") %in% names(weights))) {
    stop(
      "weights must be a list with 'center_weight' and 'neighbor_weights'",
      call. = FALSE
    )
  }
  if (!is.list(neighbors) || length(neighbors) == 0 ||
      !all(vapply(neighbors, is.list, logical(1)))) {
    stop(
      "neighbors must be a non-empty list of per-order neighbour lists",
      call. = FALSE
    )
  }
  if (length(weights$neighbor_weights) != length(neighbors)) {
    stop(
      "length(neighbor_weights) (", length(weights$neighbor_weights), ") ",
      "must equal length(neighbors) (", length(neighbors), ")",
      call. = FALSE
    )
  }
  if (length(var_names) != length(variable_values)) {
    stop("var_names length must match variable_values length", call. = FALSE)
  }
  invisible(TRUE)
}

#' Ordinal suffix for an integer (internal)
#' @keywords internal
ordinal_suffix <- function(n) {
  if (n %% 100 %in% 11:13) return("th")
  switch(as.character(n %% 10),
         "1" = "st",
         "2" = "nd",
         "3" = "rd",
         "th")
}

#' C++ wrapper for 2-order smoothing (internal)
#' @keywords internal
process_district_all_vars_wrapper <- function(variable_values, first_neighbors,
                                              second_neighbors, weights,
                                              hex_indices, var_names) {
  tryCatch(
    process_district_all_vars(
      variable_values, first_neighbors, second_neighbors,
      weights, hex_indices, var_names
    ),
    error = function(e) {
      stop("C++ function process_district_all_vars failed: ", e$message,
           call. = FALSE)
    }
  )
}

#' C++ wrapper for N-order smoothing (internal)
#' @keywords internal
process_district_all_vars_n_orders_wrapper <- function(variable_values, neighbors,
                                                       weights, hex_indices,
                                                       var_names) {
  tryCatch(
    process_district_all_vars_n_orders(
      variable_values, neighbors, weights, hex_indices, var_names
    ),
    error = function(e) {
      stop("C++ function process_district_all_vars_n_orders failed: ",
           e$message, call. = FALSE)
    }
  )
}
