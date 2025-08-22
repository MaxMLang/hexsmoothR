#' R fallback smoothing function (internal)
#' @keywords internal
#' @export
smooth_variables_r_fallback <- function(variable_values, neighbors, weights, hex_indices, var_names) {
  # Ensure we have the correct weights structure
  if (!is.list(weights) || !("center_weight" %in% names(weights)) || !("neighbor_weights" %in% names(weights))) {
    stop("Weights must be a list with 'center_weight' and 'neighbor_weights' elements")
  }
  
  # Ensure we have the correct neighbors structure
  if (!is.list(neighbors) || length(neighbors) == 0) {
    stop("Neighbors must be a non-empty list")
  }
  
  n_vars <- length(variable_values)
  n_process <- length(hex_indices)
  n_orders <- length(neighbors)
  
  result <- list()
  center_weight <- weights$center_weight
  neighbor_weights <- weights$neighbor_weights
  
  for (v in seq_len(n_vars)) {
    values <- variable_values[[v]]
    
    # Initialize result arrays for each order
    order_means <- vector("list", n_orders)
    for (order in 1:n_orders) {
      order_means[[order]] <- rep(NA_real_, n_process)
    }
    
    # Initialize weighted results
    weighted_combined <- rep(NA_real_, n_process)
    
    for (i in seq_len(n_process)) {
      hex_idx <- hex_indices[i]
      # Defensive: check bounds
      if (hex_idx < 1 || hex_idx > length(neighbors[[1]])) next
      
      # Calculate means for each order
      for (order in 1:n_orders) {
        neighs <- neighbors[[order]][[hex_idx]]
        if (length(neighs) > 0) {
          order_vals <- values[neighs]
          order_means[[order]][i] <- mean(order_vals, na.rm = TRUE)
        }
      }
      
      # Calculate weighted combined result
      weighted_sum <- 0.0
      weight_sum <- 0.0
      
             # Add center cell
       if (!is.na(values[hex_idx])) {
         weighted_sum <- weighted_sum + values[hex_idx] * center_weight
         weight_sum <- weight_sum + center_weight
       }
      
             # Add neighbors from all orders
       for (order in 1:n_orders) {
         neighs <- neighbors[[order]][[hex_idx]]
         if (length(neighs) > 0) {
           for (neigh_idx in neighs) {
             if (!is.na(values[neigh_idx])) {
               weighted_sum <- weighted_sum + values[neigh_idx] * neighbor_weights[order]
               weight_sum <- weight_sum + neighbor_weights[order]
             }
           }
         }
       }
      
      if (weight_sum > 0) {
        weighted_combined[i] <- weighted_sum / weight_sum
      }
    }
    
    # Create result structure
    var_result <- list(
      raw = weighted_combined,  # Combined result for backward compatibility
      weighted_combined = weighted_combined
    )
    
    # Add order-specific results
    for (order in 1:n_orders) {
      ordinal_suffix <- if (order == 1) "st" else if (order == 2) "nd" else if (order == 3) "rd" else "th"
      var_result[[paste0("neighbors_", order, ordinal_suffix)]] <- order_means[[order]]
    }
    

    
    result[[var_names[v]]] <- var_result
  }
  result
}

#' Apply spatial smoothing to variables using hexagonal grid topology
#' 
#' This function applies spatial smoothing to variables using a hexagonal grid topology.
#' It can be used independently after data extraction, or as part of the complete workflow.
#' The function automatically uses the C++ implementation when available and falls back to R if needed.
#' 
#' @param variable_values List of numeric vectors containing variable values for each grid cell
#' @param neighbors List of neighbor indices for each order (from compute_topology)
#' @param weights List containing center_weight and neighbor_weights for each order
#' @param hex_indices Vector of hexagon indices to process (default: all cells)
#' @param var_names Character vector of variable names
#' @return List containing smoothing results for each variable with the following components:
#'   \itemize{
#'     \item raw: Weighted average of center cell and all neighbors
#'     \item neighbors_Nst: Mean of neighbors at order N
#'     \item weighted_combined: Weighted average of center cell and all neighbors
#'   }
#' @export
#' @examples
#' # After creating a grid and computing topology
#' # grid_sf <- create_grid(study_area, cell_size = 10000)
#' # topology <- compute_topology(grid_sf, neighbor_orders = 3)  # 3 orders
#' 
#' # Extract your variables (example)
#' # variable_values <- list(
#' #   ndvi = c(0.5, 0.6, 0.4, 0.7, 0.3),
#' #   elevation = c(100, 120, 90, 140, 80)
#' # )
#' 
#' # Apply smoothing with N orders
#' # smoothing_results <- smooth_variables(
#' #   variable_values = variable_values,
#' #   neighbors = topology$neighbors,
#' #   weights = topology$weights,
#' #   var_names = c("ndvi", "elevation")
#' # )
smooth_variables <- function(variable_values, neighbors, weights, 
                           hex_indices = NULL, var_names) {
  
  # Validate required parameters
  if (is.null(neighbors) || is.null(weights)) {
    stop("Both 'neighbors' and 'weights' parameters must be provided")
  }
  
  # Handle case where user might have passed the entire topology object instead of just weights
  if (is.list(weights) && "weights" %in% names(weights)) {
    # User passed topology object, extract the weights part
    weights <- weights$weights
    cat("Note: Extracted weights from topology object\n")
  }
  
  # Handle case where user might have passed the entire topology object instead of just neighbors
  if (is.list(neighbors) && "neighbors" %in% names(neighbors)) {
    # User passed topology object, extract the neighbors part
    neighbors <- neighbors$neighbors
    cat("Note: Extracted neighbors from topology object\n")
  }
  
  # Validate weights structure
  if (!is.list(weights) || !("center_weight" %in% names(weights)) || !("neighbor_weights" %in% names(weights))) {
    stop("Weights must be a list with 'center_weight' and 'neighbor_weights' elements. Received: ", utils::str(weights))
  }
  
  # Validate neighbors structure
  if (!is.list(neighbors) || length(neighbors) == 0) {
    stop("Neighbors must be a non-empty list. Received: ", utils::str(neighbors))
  }
  
  # Set default hex_indices if not provided
  if (is.null(hex_indices)) {
    hex_indices <- 1:length(neighbors[[1]])
  }
  
  # Try to use C++ implementation if available
  result <- tryCatch({
    cat("Using C++ implementation for N-order smoothing\n")
    process_district_all_vars_n_orders(variable_values, neighbors, weights, hex_indices, var_names)
  }, error = function(e) {
    # Use R fallback if C++ fails
    cat("C++ implementation failed, falling back to R implementation: ", e$message, "\n")
    smooth_variables_r_fallback(variable_values, neighbors, weights, hex_indices, var_names)
  })
  
  return(result)
}

#' C++ wrapper for 2-order smoothing (internal)
#' @keywords internal
#' @export
process_district_all_vars_wrapper <- function(variable_values, first_neighbors, second_neighbors, weights, hex_indices, var_names) {
  tryCatch({
    return(process_district_all_vars(variable_values, first_neighbors, second_neighbors, weights, hex_indices, var_names))
  }, error = function(e) {
    stop("C++ function process_district_all_vars failed: ", e$message)
  })
}

#' C++ wrapper for N-order smoothing (internal)
#' @keywords internal
#' @export
process_district_all_vars_n_orders_wrapper <- function(variable_values, neighbors, weights, hex_indices, var_names) {
  tryCatch({
    return(process_district_all_vars_n_orders(variable_values, neighbors, weights, hex_indices, var_names))
  }, error = function(e) {
    stop("C++ function process_district_all_vars_n_orders failed: ", e$message)
  })
} 