#' Get appropriate UTM CRS for a study area
#' 
#' Automatically determines the appropriate UTM coordinate reference system (CRS) 
#' for a given study area based on its centroid location.
#' 
#' @param study_area sf object representing the study area
#' @return Character string with the UTM CRS (e.g., "EPSG:32630" for UTM zone 30N)
#' @importFrom sf st_centroid st_geometry st_coordinates
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' 
#' # Create a simple study area
#' study_area <- st_sf(geometry = st_sfc(
#'   st_polygon(list(matrix(c(-5, 35, 5, 35, 5, 45, -5, 45, -5, 35), ncol = 2, byrow = TRUE))),
#'   crs = 4326
#' ))
#' 
#' # Get appropriate UTM CRS
#' utm_crs <- get_utm_crs(study_area)
#' print(utm_crs)  # Should return something like "EPSG:32630"
#' }
get_utm_crs <- function(study_area) {
  if (!inherits(study_area, "sf")) {
    stop("study_area must be an sf object")
  }
  
  # Get centroid
  centroid <- st_centroid(st_geometry(study_area))
  coords <- st_coordinates(centroid)
  lon <- coords[1, 1]
  lat <- coords[1, 2]
  
  # Validate coordinates
  if (lon < -180 || lon > 180) {
    stop("Longitude must be between -180 and 180 degrees")
  }
  if (lat < -90 || lat > 90) {
    stop("Latitude must be between -90 and 90 degrees")
  }
  
  # Calculate UTM zone
  utm_zone <- floor((lon + 180) / 6) + 1
  
  # Handle edge cases
  if (lon >= 180) utm_zone <- 60
  if (lon < -180) utm_zone <- 1
  
  # Determine hemisphere
  if (lat >= 0) {
    utm_crs <- paste0("EPSG:", 32600 + utm_zone)
  } else {
    utm_crs <- paste0("EPSG:", 32700 + utm_zone)
  }
  
  return(utm_crs)
}

#' Convert between hexagon measurements
#' 
#' Helper functions to convert between different hexagon measurements:
#' flat-to-flat distance, edge length, and circumradius (center to vertex).
#' 
#' Convert between hexagon measurements
#' 
#' These functions convert between different hexagon measurements:
#' - `hex_flat_to_edge`: Convert flat-to-flat distance to edge length
#' - `hex_flat_to_circumradius`: Convert flat-to-flat distance to circumradius
#' - `hex_edge_to_flat`: Convert edge length to flat-to-flat distance  
#' - `hex_circumradius_to_flat`: Convert circumradius to flat-to-flat distance
#' 
#' @param flat_to_flat Flat-to-flat distance (distance between opposite edges)
#' @param edge_length Edge length of the hexagon
#' @param circumradius Circumradius of the hexagon
#' @return Numeric value in the same units as input
#' @export
#' @examples
#' # Convert 1000m flat-to-flat to edge length
#' hex_flat_to_edge(1000)  # Returns ~577.35m
#' 
#' # Convert 1000m flat-to-flat to circumradius  
#' hex_flat_to_circumradius(1000)  # Returns ~577.35m
#' 
#' # Convert 577.35m edge length to flat-to-flat
#' hex_edge_to_flat(577.35)  # Returns ~1000m
#' 
#' # Convert 577.35m circumradius to flat-to-flat
#' hex_circumradius_to_flat(577.35)  # Returns ~1000m
hex_flat_to_edge <- function(flat_to_flat) {
  flat_to_flat / sqrt(3)
}

#' @rdname hex_flat_to_edge
#' @export
hex_flat_to_circumradius <- function(flat_to_flat) {
  flat_to_flat / sqrt(3)
}

#' @rdname hex_flat_to_edge  
#' @export
hex_edge_to_flat <- function(edge_length) {
  edge_length * sqrt(3)
}

#' @rdname hex_flat_to_edge
#' @export
hex_circumradius_to_flat <- function(circumradius) {
  circumradius * sqrt(3)
}

#' Create hexagonal or square grids for spatial analysis
#'
#' Creates regular hexagonal or square grids over a study area. The function automatically
#' handles coordinate system transformations and ensures proper grid alignment.
#'
#' @details 
#' **CRITICAL: Cell Size Units**
#' 
#' The `cell_size` parameter units depend on the coordinate reference system (CRS):
#' 
#' - **Projected CRS (UTM, State Plane, etc.)**: `cell_size` in **METERS**
#'   - Example: `cell_size = 20000` creates 20km hexagons
#'   - Example: `cell_size = 1000` creates 1km hexagons
#' 
#' - **Geographic CRS (WGS84, NAD83, etc.)**: `cell_size` in **DEGREES**
#'   - Example: `cell_size = 0.1` creates ~11km hexagons
#'   - Example: `cell_size = 0.01` creates ~1.1km hexagons
#' 
#' **Recommended Practice:** Always use projected CRS (UTM) for real-world analysis
#' to ensure accurate area calculations and distance measurements.
#' 
#' **Hexagon Measurements:**
#' - `cell_size` represents the **flat-to-flat distance** (width between opposite edges)
#' - Creates **pointy-topped hexagons** (flat sides horizontal)
#' - This is the most intuitive measure as it represents the "width" of each hexagon
#' - Use helper functions `hex_flat_to_edge()` etc. to convert between measurements
#'
#' @param study_area sf object containing polygonal geometries defining the study area.
#'   Can be in any CRS - will be transformed to `projection_crs` for grid creation.
#' @param cell_size Grid cell size. **UNITS DEPEND ON CRS:**
#'   - **Projected CRS (recommended)**: size in METERS
#'     - For hexagons: flat-to-flat distance (width between opposite edges)
#'     - For squares: side length
#'   - **Geographic CRS**: size in DEGREES
#'     - Use with caution for large areas due to distortion
#' @param type Character string: "hexagonal" (default) or "square"
#' @param projection_crs CRS for grid creation (default UTM 36N). Should be a projected
#'   coordinate system for accurate measurements. Use `get_utm_crs(study_area)` to
#'   automatically determine appropriate UTM zone.
#' @param id_column Optional column name in `study_area` for creating separate grids
#'   for each unique value. If provided, returns a named list of grids.
#' @param return_crs CRS for output grid (default WGS84). Grid will be transformed
#'   to this CRS after creation.
#' @param check_size Whether to warn if grid will be very large (default TRUE)
#' @param max_cells Maximum number of cells allowed before stopping (default 1000000).
#'   Set to NULL to disable this safety check entirely.
#' @return sf object containing the grid, or named list of sf objects if `id_column` is used.
#'   Each grid cell contains:
#'   \describe{
#'     \item{grid_id}{Unique identifier for each cell}
#'     \item{grid_index}{Sequential index (1 to number of cells)}
#'     \item{geometry}{Polygon geometry of the cell}
#'   }
#' @importFrom sf st_geometry_type st_make_grid st_transform st_bbox st_crop st_intersection
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' 
#' # Create a simple study area
#' study_area <- st_sf(geometry = st_sfc(
#'   st_polygon(list(matrix(c(-5, 35, 5, 35, 5, 45, -5, 45, -5, 35), ncol = 2, byrow = TRUE))),
#'   crs = 4326
#' ))
#' 
#' # Create 1km hexagonal grid (flat-to-flat distance = 1000m)
#' hex_grid <- create_grid(study_area, cell_size = 1000, type = "hexagonal")
#' 
#' # Create 2km square grid (side length = 2000m)  
#' square_grid <- create_grid(study_area, cell_size = 2000, type = "square")
#' 
#' # Create high-resolution grid with custom safety limit
#' high_res_grid <- create_grid(study_area, cell_size = 100, type = "hexagonal", max_cells = 2000000)
#' }
create_grid <- function(study_area, 
                       cell_size, 
                       type = c("hexagonal", "square"), 
                       projection_crs = 32636,
                       id_column = NULL,
                       return_crs = 4326,
                       check_size = TRUE,
                       max_cells = 1000000) {
  # quick checks
  if (!inherits(study_area, "sf")) stop("study_area must be sf")
  if (cell_size <= 0) stop("cell_size must be positive")
  type <- match.arg(type)
  grid_size <- cell_size
  # handle multiple areas
  if (!is.null(id_column)) {
    if (!id_column %in% names(study_area)) stop("id_column not found")
    unique_areas <- unique(study_area[[id_column]])
    grids <- list()
    for (area_id in unique_areas) {
      cat("Processing area:", area_id, "\n")
      area_sf <- study_area[study_area[[id_column]] == area_id, ]
      grids[[as.character(area_id)]] <- create_single_grid(
        area_sf, grid_size, type, projection_crs, return_crs, area_id, check_size, max_cells
      )
    }
    return(grids)
  } else {
    return(create_single_grid(study_area, grid_size, type, projection_crs, return_crs, NULL, check_size, max_cells))
  }
}

#' Find optimal hexagonal cell size for target number of cells
#' 
#' Uses binary search to find the hexagonal cell size (flat-to-flat distance) 
#' that produces approximately the target number of hexagons.
#' 
#' @param study_area sf object (polygons)
#' @param target_cells target number of hexagons  
#' @param cell_size_min minimum cell size to try (in map units)
#' @param cell_size_max maximum cell size to try (in map units)
#' @param tol tolerance for convergence (fraction of target_cells)
#' @param max_iter maximum iterations for binary search
#' @return cell size (flat-to-flat distance) that gives closest to target_cells hexagons
#' @export
find_hex_cell_size_for_target_cells <- function(study_area, target_cells, cell_size_min = 0.001, cell_size_max = 10, tol = 0.05, max_iter = 20) {
  if (!inherits(study_area, "sf")) stop("study_area must be sf")
  if (target_cells <= 0) stop("target_cells must be positive")
  lower <- cell_size_min
  upper <- cell_size_max
  best_size <- NA
  best_diff <- Inf
  n_cells_min <- tryCatch({
    grid <- create_grid(study_area, cell_size = lower, type = "hexagonal", check_size = FALSE)
    nrow(grid)
  }, error = function(e) NA)
  n_cells_max <- tryCatch({
    grid <- create_grid(study_area, cell_size = upper, type = "hexagonal", check_size = FALSE)
    nrow(grid)
  }, error = function(e) NA)
  cat(sprintf("[find_hex_cell_size_for_target_cells] n_cells at min (%.4f): %s\n", lower, n_cells_min))
  cat(sprintf("[find_hex_cell_size_for_target_cells] n_cells at max (%.4f): %s\n", upper, n_cells_max))
  for (i in 1:max_iter) {
    mid <- (lower + upper) / 2
    n_cells <- tryCatch({
      grid <- create_grid(study_area, cell_size = mid, type = "hexagonal", check_size = FALSE)
      nrow(grid)
    }, error = function(e) NA)
    if (is.na(n_cells)) {
      lower <- mid
      next
    }
    diff <- abs(n_cells - target_cells)
    if (diff < best_diff) {
      best_diff <- diff
      best_size <- mid
    }
    if (diff / target_cells < tol) return(mid)
    if (n_cells > target_cells) {
      lower <- mid
    } else {
      upper <- mid
    }
  }
  if (is.na(best_size)) {
    cat("[find_hex_cell_size_for_target_cells] Could not find a valid cell size. Returning max cell size.\n")
    return(cell_size_max)
  }
  cat(sprintf("[find_hex_cell_size_for_target_cells] Closest found: %.4f (diff: %d)\n", best_size, best_diff))
  return(best_size)
}

#' Actually make the grid for one area
#' @keywords internal
create_single_grid <- function(study_area, grid_size, type, projection_crs, return_crs, area_id = NULL, check_size = TRUE, max_cells = 1000000) {
  # checks
  if (!inherits(study_area, "sf")) stop("study_area must be sf")
  if (nrow(study_area) == 0) stop("study_area is empty")
  if (grid_size <= 0) stop("grid_size must be positive")
  geom_types <- sf::st_geometry_type(study_area)
  if (!all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))) stop("study_area must be polygons")
  # project for grid creation
  study_area_proj <- tryCatch({
    sf::st_transform(study_area, projection_crs)
  }, error = function(e) stop("Failed to transform study_area: ", e$message))
  bbox <- tryCatch({
    sf::st_bbox(study_area_proj)
  }, error = function(e) stop("Failed to get bbox: ", e$message))
  width <- bbox[3] - bbox[1]
  height <- bbox[4] - bbox[2]
  if (width <= 0 || height <= 0) stop("Invalid extent")
  # just try to make the grid
  cat("  Creating", type, "grid...\n")
  grid <- tryCatch({
    sf::st_sf(sf::st_make_grid(study_area_proj, cellsize = grid_size, square = (type == "square")))
  }, error = function(e) {
    if (grepl("memory", e$message, ignore.case = TRUE)) stop("Memory error: use bigger cell_size")
    else stop("Error making grid: ", e$message)
  })
  if (check_size && nrow(grid) > 100000) warning("Grid is very large (", nrow(grid), ")")
  if (check_size && !is.null(max_cells) && nrow(grid) > max_cells) stop("Grid is too large (", nrow(grid), "). Increase max_cells parameter or use larger cell_size.")
  # crop and clip
  cat("  Cropping grid...\n")
  grid_cropped <- tryCatch({
    sf::st_crop(grid, study_area_proj)
  }, error = function(e) stop("Error cropping: ", e$message))
  rm(grid); gc(FALSE)
  cat("  Clipping grid...\n")
  grid_clipped <- tryCatch({
    sf::st_intersection(grid_cropped, study_area_proj)
  }, error = function(e) stop("Error clipping: ", e$message))
  rm(grid_cropped); gc(FALSE)
  # back to output CRS
  cat("  Transforming CRS...\n")
  grid_final <- tryCatch({
    sf::st_transform(grid_clipped, crs = return_crs)
  }, error = function(e) stop("Error transforming: ", e$message))
  rm(grid_clipped); gc(FALSE)
  # only polygons
  grid_final <- grid_final[sf::st_geometry_type(grid_final) %in% c("POLYGON", "MULTIPOLYGON"), ]
  if (nrow(grid_final) == 0) stop("No valid grid cells created")
  # add IDs
  if (!is.null(area_id)) {
    grid_final$area_id <- area_id
    grid_final$grid_id <- paste0(area_id, "_", seq_len(nrow(grid_final)))
  } else {
    grid_final$grid_id <- paste0("grid_", seq_len(nrow(grid_final)))
  }
  grid_final$grid_index <- seq_len(nrow(grid_final))
  cat("  Grid created with", nrow(grid_final), "cells\n")
  return(grid_final)
} 