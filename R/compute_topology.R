#' Compute spatial topology and weights for hexagonal smoothing
#'
#' Finds spatial neighbours and computes Gaussian-based weights for smoothing.
#' Neighbours are computed up to `neighbor_orders` orders (1st-order = touching,
#' 2nd-order = neighbours of neighbours, etc.). Weights decay with order using a
#' Gaussian kernel.
#'
#' @section Weight computation:
#' - `sigma = avg_distance * adaptive_sigma_factor` (auto-bandwidth) when
#'   `sigma` is `NULL`.
#' - For order N: weight ~ `exp(-N^2 / (2 * sigma^2))`.
#' - All weights (including `center_weight`) are normalised to sum to 1.
#'
#' @param grid sf object with polygonal geometries (or a list of sf objects).
#' @param projection_crs CRS used for distance calculations. Default `NULL`
#'   selects an appropriate UTM zone via [`get_utm_crs()`].
#' @param neighbor_orders Number of neighbour orders (positive integer).
#' @param sigma Gaussian bandwidth. `NULL` (default) auto-computes from
#'   `avg_distance` and `adaptive_sigma_factor`.
#' @param center_weight Weight for the centre cell.
#' @param neighbor_weights_param Optional list of length `neighbor_orders` with
#'   per-order weights (overrides Gaussian computation).
#' @param adaptive_sigma_factor Scaling factor for auto-sigma.
#' @param sample_size Cells to sample for the average-distance estimate.
#' @return List (or named list of lists if `grid` was a list) containing
#'   `neighbors`, `weights`, `avg_distance`, `sigma`, `grid_ids`,
#'   `grid_indices`, `neighbor_orders`.
#' @importFrom sf st_transform st_touches st_centroid st_distance st_crs
#' @export
compute_topology <- function(grid,
                             projection_crs = NULL,
                             neighbor_orders = 2,
                             sigma = NULL,
                             center_weight = 1.0,
                             neighbor_weights_param = NULL,
                             adaptive_sigma_factor = 0.5,
                             sample_size = 100) {
  if (!is.numeric(neighbor_orders) ||
      length(neighbor_orders) != 1L ||
      neighbor_orders < 1 ||
      neighbor_orders != round(neighbor_orders)) {
    stop("neighbor_orders must be a positive integer", call. = FALSE)
  }
  neighbor_orders <- as.integer(neighbor_orders)

  if (!is.null(neighbor_weights_param)) {
    if (!is.list(neighbor_weights_param) ||
        length(neighbor_weights_param) != neighbor_orders) {
      stop(
        "neighbor_weights_param must be a list with length equal to neighbor_orders",
        call. = FALSE
      )
    }
    if (!all(vapply(neighbor_weights_param, is.numeric, logical(1)))) {
      stop("All elements in neighbor_weights_param must be numeric", call. = FALSE)
    }
    neighbor_weights_param <- unlist(neighbor_weights_param)
  }

  if (inherits(grid, "sf")) {
    grid_list <- list(grid = grid)
    single_grid <- TRUE
  } else if (is.list(grid) && all(vapply(grid, inherits, logical(1), "sf"))) {
    grid_list <- grid
    if (is.null(names(grid_list))) {
      names(grid_list) <- paste0("grid_", seq_along(grid_list))
    }
    single_grid <- FALSE
  } else {
    stop("grid must be an sf object or a list of sf objects", call. = FALSE)
  }

  topology <- vector("list", length(grid_list))
  names(topology) <- names(grid_list)

  for (area_name in names(grid_list)) {
    grid_sf <- grid_list[[area_name]]
    hex_msg(
      "Computing topology for area: ", area_name, " (",
      nrow(grid_sf), " cells, ", neighbor_orders, " orders)\n"
    )

    if (nrow(grid_sf) == 0) {
      topology[area_name] <- list(NULL)
      next
    }

    if (!"grid_id" %in% names(grid_sf)) {
      grid_sf$grid_id <- paste0(area_name, "_", seq_len(nrow(grid_sf)))
    }
    if (!"grid_index" %in% names(grid_sf)) {
      grid_sf$grid_index <- seq_len(nrow(grid_sf))
    }

    proj_crs <- projection_crs
    if (is.null(proj_crs)) {
      proj_crs <- tryCatch(get_utm_crs(grid_sf), error = function(e) sf::st_crs(grid_sf))
    }
    grid_proj <- tryCatch(
      sf::st_transform(grid_sf, proj_crs),
      error = function(e) {
        warning("Could not transform grid for ", area_name, ": ", e$message, call. = FALSE)
        grid_sf
      }
    )

    neighbors <- compute_neighbors_n_orders(grid_proj, neighbor_orders)
    if (is.null(neighbors)) {
      hex_msg("  Warning: no neighbours found\n")
      topology[area_name] <- list(NULL)
      next
    }

    avg_dist <- estimate_avg_neighbor_distance(grid_proj, neighbors[[1]], sample_size)
    if (is.na(avg_dist)) {
      hex_msg("  Warning: could not calculate distances\n")
      topology[area_name] <- list(NULL)
      next
    }

    weights <- compute_weights_n_orders(
      neighbor_orders, sigma, center_weight, neighbor_weights_param,
      avg_dist, adaptive_sigma_factor
    )

    topology[[area_name]] <- list(
      neighbors = neighbors,
      avg_distance = avg_dist,
      sigma = weights$sigma,
      weights = weights,
      grid_ids = grid_sf$grid_id,
      grid_indices = grid_sf$grid_index,
      neighbor_orders = neighbor_orders
    )

    hex_msg(
      "  Computed: avg_dist = ", round(avg_dist, 2),
      ", orders = ", neighbor_orders, "\n"
    )
  }

  if (single_grid && length(topology) == 1) return(topology[[1]])
  topology
}

#' Estimate average distance between first-order neighbours (internal)
#' @keywords internal
estimate_avg_neighbor_distance <- function(grid_proj, first_neighbors, sample_size) {
  n_cells <- nrow(grid_proj)
  if (n_cells < 2) return(NA_real_)

  has_neigh <- vapply(first_neighbors, function(x) length(x) > 0, logical(1))
  candidates <- which(has_neigh)
  if (length(candidates) == 0) return(NA_real_)

  k <- min(sample_size, length(candidates))
  sample_indices <- if (length(candidates) == 1) candidates else sample(candidates, k)

  centroids <- suppressWarnings(sf::st_centroid(grid_proj))
  distances <- numeric(0)

  for (i in sample_indices) {
    nb <- first_neighbors[[i]]
    if (length(nb) > 3) nb <- sample(nb, 3)
    if (length(nb) == 0) next
    d <- as.numeric(sf::st_distance(centroids[i, ], centroids[nb, ]))
    distances <- c(distances, d)
  }

  if (length(distances) == 0) return(NA_real_)
  mean(distances, na.rm = TRUE)
}

#' Compute neighbours for N orders (internal, vectorised set-based BFS)
#' @keywords internal
compute_neighbors_n_orders <- function(grid_proj, neighbor_orders) {
  n_cells <- nrow(grid_proj)

  first_neighbors <- tryCatch(
    sf::st_touches(grid_proj),
    error = function(e) {
      hex_msg("  Error computing first-order neighbours: ", e$message, "\n")
      NULL
    }
  )
  if (is.null(first_neighbors) || length(first_neighbors) == 0) return(NULL)

  neighbors <- vector("list", neighbor_orders)
  neighbors[[1]] <- lapply(first_neighbors, as.integer)

  if (neighbor_orders >= 2) {
    self_seq <- seq_len(n_cells)
    seen <- vector("list", n_cells)
    for (i in self_seq) {
      seen[[i]] <- as.integer(c(i, neighbors[[1]][[i]]))
    }
    prev_frontier <- neighbors[[1]]

    for (order in 2:neighbor_orders) {
      next_frontier <- vector("list", n_cells)
      for (i in self_seq) {
        front <- prev_frontier[[i]]
        if (length(front) == 0) {
          next_frontier[[i]] <- integer(0)
          next
        }
        cand <- unlist(neighbors[[1]][front], use.names = FALSE)
        if (length(cand) == 0) {
          next_frontier[[i]] <- integer(0)
          next
        }
        cand <- unique.default(cand)
        new_cells <- cand[!cand %in% seen[[i]]]
        next_frontier[[i]] <- new_cells
        if (length(new_cells)) {
          seen[[i]] <- c(seen[[i]], new_cells)
        }
      }
      neighbors[[order]] <- next_frontier
      prev_frontier <- next_frontier
    }
  }

  neighbors
}

#' Compute weights for N orders (internal)
#' @keywords internal
compute_weights_n_orders <- function(neighbor_orders, sigma, center_weight,
                                     weights_input, avg_dist,
                                     adaptive_sigma_factor) {
  if (is.null(sigma)) sigma <- avg_dist * adaptive_sigma_factor

  if (is.null(weights_input)) {
    orders <- seq_len(neighbor_orders)
    neighbor_weights <- exp(-orders^2 / (2 * sigma^2))
  } else {
    if (is.list(weights_input)) {
      neighbor_weights <- unlist(weights_input)
    } else {
      neighbor_weights <- weights_input
    }
    if (!is.numeric(neighbor_weights)) {
      stop("weights_input must be numeric or a list of numerics", call. = FALSE)
    }
  }

  total_weight <- center_weight + sum(neighbor_weights)
  if (total_weight > 0) {
    center_weight <- center_weight / total_weight
    neighbor_weights <- neighbor_weights / total_weight
  }

  list(
    center_weight = center_weight,
    neighbor_weights = neighbor_weights,
    sigma = sigma
  )
}
