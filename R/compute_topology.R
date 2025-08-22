#' Compute spatial topology and weights for hexagonal smoothing
#'
#' Finds spatial neighbors and computes Gaussian-based weights for smoothing. 
#' This function determines which grid cells are neighbors at specified orders
#' (e.g., first-order neighbors are touching, second-order are neighbors of neighbors),
#' then calculates appropriate weights based on spatial distance using Gaussian decay.
#'
#' @details 
#' **Weight Computation:**
#' 
#' The function uses Gaussian weights based on spatial distance:
#' - `sigma = avg_distance * adaptive_sigma_factor` (adaptive bandwidth)
#' - Weights decay exponentially with neighbor order: `exp(-order^2 / (2 * sigma^2))`
#' - All weights are normalized so they sum to 1
#' 
#' **How weights are applied in smoothing:**
#' - Each individual neighbor gets the full weight (not averaged)
#' - Center cell: `value * center_weight`
#' - Each neighbor at order N: `value * neighbor_weight_order_N`
#' - Final value = sum of all weighted values / sum of all weights used
#' 
#' This means cells with more neighbors get more influence from their neighborhood,
#' which is realistic for spatial analysis.
#'
#' @param grid sf object with polygonal geometries (hexagons or squares) or list of sf objects.
#'   Must contain columns `grid_id` and `grid_index`, or they will be created automatically.
#' @param projection_crs CRS for distance calculations (default UTM 36N). Should be a projected
#'   coordinate system (like UTM) for accurate distance measurements in meters.
#' @param neighbor_orders Number of neighbor orders to compute (default 2).
#'   Higher orders provide more extensive smoothing but increase computation time and memory usage.
#' @param sigma Gaussian bandwidth parameter in meters (NULL = auto-computed as avg_distance * adaptive_sigma_factor).
#'   Controls how quickly weights decay with distance.
#' @param center_weight Weight for the center cell (default 1.0). Higher values give more 
#'   influence to the original cell value.
#' @param neighbor_weights_param Optional list of weights for each neighbor order (overrides Gaussian computation).
#'   If specified, must have length equal to `neighbor_orders`. Example: `list(0.5, 0.25)` for 2 orders.
#' @param adaptive_sigma_factor Scaling factor for automatic sigma computation (default 0.5).
#'   Smaller values = less smoothing, larger values = more smoothing.
#' @param sample_size Number of cells to sample for distance calculations (default 100).
#'   Used to estimate average neighbor distance efficiently.
#' @return List containing:
#' \describe{
#'   \item{neighbors}{List of lists, each containing neighbor indices for each order}
#'   \item{weights}{List with center_weight and neighbor weights for each order (all normalized)}
#'   \item{avg_distance}{Average distance between neighboring cells in meters}
#'   \item{sigma}{Gaussian bandwidth parameter used}
#'   \item{grid_ids}{Character vector of grid cell identifiers}
#'   \item{grid_indices}{Integer vector of grid cell indices}
#'   \item{neighbor_orders}{Number of neighbor orders computed}
#' }
#' @importFrom sf st_transform st_touches st_centroid st_distance
#' @export
compute_topology <- function(grid,
                           projection_crs = 32636,
                           neighbor_orders = 2,
                           sigma = NULL,
                           center_weight = 1.0,
                           neighbor_weights_param = NULL,
                           adaptive_sigma_factor = 0.5,
                           sample_size = 100) {
  # Validate neighbor_orders parameter
  if (!is.numeric(neighbor_orders) || neighbor_orders < 1 || neighbor_orders != round(neighbor_orders)) {
    stop("neighbor_orders must be a positive integer")
  }
  
  # Validate neighbor_weights_param if provided
  if (!is.null(neighbor_weights_param)) {
    if (!is.list(neighbor_weights_param) || length(neighbor_weights_param) != neighbor_orders) {
      stop("neighbor_weights_param must be a list with length equal to neighbor_orders")
    }
    if (!all(sapply(neighbor_weights_param, is.numeric))) {
      stop("All elements in neighbor_weights_param must be numeric")
    }
    # Convert list to numeric vector for internal processing
    neighbor_weights_param <- unlist(neighbor_weights_param)
  }
  
  # handle single or list
  if (inherits(grid, "sf")) {
    grid_list <- list("grid" = grid)
    single_grid <- TRUE
  } else if (is.list(grid) && all(sapply(grid, inherits, "sf"))) {
    grid_list <- grid
    single_grid <- FALSE
  } else {
    stop("grid must be an sf object or a list of sf objects")
  }
  
  topology <- list()
  for (area_name in names(grid_list)) {
    cat("Computing topology for area:", area_name, "with", neighbor_orders, "neighbor orders\n")
    grid_sf <- grid_list[[area_name]]
    
    # add IDs if missing
    if (!"grid_id" %in% names(grid_sf)) grid_sf$grid_id <- paste0(area_name, "_", 1:nrow(grid_sf))
    if (!"grid_index" %in% names(grid_sf)) grid_sf$grid_index <- 1:nrow(grid_sf)
    
    # project for distances
    grid_proj <- st_transform(grid_sf, projection_crs)
    
    # Compute neighbors for all orders
    neighbors <- compute_neighbors_n_orders(grid_proj, neighbor_orders)
    
    if (is.null(neighbors)) { 
      cat("  Warning: No neighbors found\n"); 
      topology[[area_name]] <- NULL; 
      next 
    }
    
    # get average distance between neighbors
    centroids <- st_centroid(grid_proj)
    actual_sample_size <- min(sample_size, nrow(grid_proj))
    distances <- numeric(0)
    
    # Sample more representative distances
    sample_indices <- sample(1:nrow(grid_proj), actual_sample_size)
    for (i in sample_indices) {
      if (length(neighbors[[1]][[i]]) > 0) {  # Check first-order neighbors
        # Sample multiple neighbors for better representation
        neighbor_indices <- neighbors[[1]][[i]]
        if (length(neighbor_indices) > 3) {
          neighbor_indices <- sample(neighbor_indices, 3)
        }
        for (neighbor_idx in neighbor_indices) {
          dist <- st_distance(centroids[i,], centroids[neighbor_idx,])
          distances <- c(distances, as.numeric(dist))
        }
      }
    }
    
    if (length(distances) == 0) { 
      cat("  Warning: Could not calculate distances\n"); 
      topology[[area_name]] <- NULL; 
      next 
    }
    
    avg_dist <- mean(distances, na.rm = TRUE)
    
    # Compute weights for all neighbor orders
    weights <- compute_weights_n_orders(
      neighbor_orders, sigma, center_weight, neighbor_weights_param, 
      avg_dist, adaptive_sigma_factor
    )
    
    # Create result structure
    result <- list(
      neighbors = neighbors,
      avg_distance = avg_dist,
      sigma = weights$sigma,
      weights = weights,
      grid_ids = grid_sf$grid_id,
      grid_indices = grid_sf$grid_index,
      neighbor_orders = neighbor_orders
    )
    

    
    topology[[area_name]] <- result
    cat("  Computed: avg_dist =", round(avg_dist, 2), "m, orders =", neighbor_orders, "\n")
  }
  
  if (single_grid && length(topology) == 1) return(topology[[1]])
  return(topology)
}

#' Compute neighbors for N orders (internal function)
#' @keywords internal
compute_neighbors_n_orders <- function(grid_proj, neighbor_orders) {
  n_cells <- nrow(grid_proj)
  
  # Initialize neighbors list
  neighbors <- vector("list", neighbor_orders)
  for (i in 1:neighbor_orders) {
    neighbors[[i]] <- vector("list", n_cells)
  }
  
  # Compute first-order neighbors (touching cells)
  first_neighbors <- tryCatch({ 
    st_touches(grid_proj) 
  }, error = function(e) { 
    cat("  Error computing first-order neighbors: ", e$message, "\n"); 
    return(list()) 
  })
  
  if (length(first_neighbors) == 0) return(NULL)
  
  # Store first-order neighbors
  neighbors[[1]] <- first_neighbors
  
  # Compute higher-order neighbors using breadth-first search
  for (order in 2:neighbor_orders) {
    for (i in 1:n_cells) {
      # Use breadth-first search to find cells at exactly 'order' steps away
      current_neighbors <- integer(0)
      visited <- integer(0)
      queue <- list(list(cell = i, distance = 0))
      
      while (length(queue) > 0) {
        current <- queue[[1]]
        queue <- queue[-1]
        
        if (current$distance == order) {
          # We've reached the target distance
          if (current$cell != i) {  # Don't include self
            current_neighbors <- c(current_neighbors, current$cell)
          }
        } else if (current$distance < order) {
          # Continue searching
          if (!(current$cell %in% visited)) {
            visited <- c(visited, current$cell)
            
            # Get neighbors of current cell
            if (current$cell <= length(first_neighbors) && length(first_neighbors[[current$cell]]) > 0) {
              for (neighbor in first_neighbors[[current$cell]]) {
                if (!(neighbor %in% visited)) {
                  queue <- c(queue, list(list(cell = neighbor, distance = current$distance + 1)))
                }
              }
            }
          }
        }
      }
      
      # Remove duplicates and cells from previous orders
      all_prev_orders <- integer(0)
      for (prev_order in 1:(order-1)) {
        all_prev_orders <- c(all_prev_orders, neighbors[[prev_order]][[i]])
      }
      current_neighbors <- setdiff(unique(current_neighbors), c(i, all_prev_orders))
      
      neighbors[[order]][[i]] <- current_neighbors
    }
  }
  
  return(neighbors)
}

#' Compute weights for N orders (internal function)
#' @keywords internal
compute_weights_n_orders <- function(neighbor_orders, sigma, center_weight, weights_input, 
                                   avg_dist, adaptive_sigma_factor) {
  # Set up sigma if not provided
  if (is.null(sigma)) sigma <- avg_dist * adaptive_sigma_factor
  
  # Compute weights
  if (is.null(weights_input)) {
    # Use Gaussian weights based on order
    neighbor_weights <- numeric(neighbor_orders)
    for (order in 1:neighbor_orders) {
      neighbor_weights[order] <- exp(-order^2 / (2 * sigma^2))
    }
  } else {
    # Ensure weights_input is a numeric vector
    if (is.list(weights_input)) {
      neighbor_weights <- unlist(weights_input)
    } else {
      neighbor_weights <- weights_input
    }
    if (!is.numeric(neighbor_weights)) {
      stop("weights_input must be numeric or a list of numeric values")
    }
  }
  
  # Normalize weights to ensure they sum to 1
  total_weight <- center_weight + sum(neighbor_weights)
  if (total_weight > 0) {
    center_weight <- center_weight / total_weight
    neighbor_weights <- neighbor_weights / total_weight
  }
  
  return(list(
    center_weight = center_weight,
    neighbor_weights = neighbor_weights,
    sigma = sigma
  ))
}

 