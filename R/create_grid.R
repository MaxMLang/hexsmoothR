#' Get an appropriate UTM CRS for a study area
#'
#' Determines the appropriate UTM coordinate reference system (CRS) for a given
#' study area based on the longitude/latitude of its centroid.
#'
#' @param study_area sf object representing the study area. May be in any CRS;
#'   it will be transformed to WGS84 (EPSG:4326) internally for the centroid
#'   calculation.
#' @return Character string with the UTM CRS, e.g. `"EPSG:32630"` for UTM 30N.
#' @importFrom sf st_centroid st_geometry st_coordinates st_transform st_crs st_is_longlat
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' study_area <- st_sf(geometry = st_sfc(
#'   st_polygon(list(matrix(c(-5, 35, 5, 35, 5, 45, -5, 45, -5, 35),
#'                          ncol = 2, byrow = TRUE))),
#'   crs = 4326
#' ))
#' get_utm_crs(study_area)  # "EPSG:32630"
#' }
get_utm_crs <- function(study_area) {
  if (!inherits(study_area, "sf")) {
    stop("study_area must be an sf object", call. = FALSE)
  }
  if (is.na(st_crs(study_area))) {
    stop("study_area has no CRS; assign one before calling get_utm_crs()", call. = FALSE)
  }

  if (!isTRUE(st_is_longlat(study_area))) {
    study_area <- st_transform(study_area, 4326)
  }

  centroid <- suppressWarnings(st_centroid(st_geometry(study_area)))
  coords <- st_coordinates(centroid)
  lon <- coords[1, 1]
  lat <- coords[1, 2]

  if (is.na(lon) || is.na(lat)) {
    stop("Could not compute centroid coordinates for study_area", call. = FALSE)
  }
  if (lon < -180 || lon > 180) {
    stop("Longitude must be between -180 and 180 degrees", call. = FALSE)
  }
  if (lat < -90 || lat > 90) {
    stop("Latitude must be between -90 and 90 degrees", call. = FALSE)
  }

  utm_zone <- floor((lon + 180) / 6) + 1
  if (utm_zone > 60) utm_zone <- 60   # handles lon == 180

  if (lat >= 0) {
    paste0("EPSG:", 32600 + utm_zone)
  } else {
    paste0("EPSG:", 32700 + utm_zone)
  }
}

#' Convert between hexagon measurements
#'
#' Helper functions to convert between different hexagon measurements.
#' For a regular hexagon, the circumradius equals the edge length, so the
#' "edge" and "circumradius" helpers are mathematically identical (provided
#' for clarity at the call site).
#'
#' - `hex_flat_to_edge()`: flat-to-flat distance to edge length
#' - `hex_flat_to_circumradius()`: flat-to-flat distance to circumradius
#' - `hex_edge_to_flat()`: edge length to flat-to-flat distance
#' - `hex_circumradius_to_flat()`: circumradius to flat-to-flat distance
#'
#' @param flat_to_flat Flat-to-flat distance (between opposite edges).
#' @param edge_length Edge length of the hexagon.
#' @param circumradius Circumradius (centre to vertex).
#' @return Numeric value in the same units as the input.
#' @export
#' @examples
#' hex_flat_to_edge(1000)          # ~577.35
#' hex_flat_to_circumradius(1000)  # ~577.35
#' hex_edge_to_flat(577.35)        # ~1000
#' hex_circumradius_to_flat(577.35)# ~1000
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
#' Creates regular hexagonal or square grids over a study area. The function
#' automatically handles coordinate-system transformations and ensures proper
#' grid alignment.
#'
#' @details
#' **Cell size units depend on the projection CRS:**
#' - Projected CRS (UTM, etc.): `cell_size` is in **metres**.
#' - Geographic CRS (WGS84): `cell_size` is in **degrees**.
#'
#' Use a projected CRS for real-world analysis. If `projection_crs` is `NULL`
#' (the default), an appropriate UTM zone is chosen automatically using
#' [`get_utm_crs()`].
#'
#' Hexagons are pointy-topped; `cell_size` is the flat-to-flat distance
#' (width between opposite edges).
#'
#' @param study_area sf object containing polygonal geometries defining the
#'   study area. May be in any CRS.
#' @param cell_size Grid cell size (see Details for units).
#' @param type `"hexagonal"` (default) or `"square"`.
#' @param projection_crs CRS used for grid construction. Default `NULL` selects
#'   an appropriate UTM zone via [`get_utm_crs()`]. Pass an EPSG code or CRS
#'   string to override.
#' @param id_column Optional column name in `study_area` for creating separate
#'   grids per unique value. If supplied, returns a named list of grids.
#' @param return_crs CRS for the output grid (default WGS84). The grid is
#'   transformed to this CRS after creation.
#' @param check_size Whether to warn when the grid is very large (default TRUE).
#' @param max_cells Maximum number of cells allowed before stopping
#'   (default 1,000,000). Set to `NULL` to disable.
#' @return sf object containing the grid (or named list of sf objects when
#'   `id_column` is given). Each grid cell carries `grid_id` and `grid_index`.
#' @importFrom sf st_geometry_type st_make_grid st_transform st_bbox st_crop
#'   st_intersection st_crs
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' study_area <- st_sf(geometry = st_sfc(
#'   st_polygon(list(matrix(c(-5, 35, 5, 35, 5, 45, -5, 45, -5, 35),
#'                          ncol = 2, byrow = TRUE))),
#'   crs = 4326
#' ))
#' hex_grid <- create_grid(study_area, cell_size = 1000, type = "hexagonal")
#' }
create_grid <- function(study_area,
                        cell_size,
                        type = c("hexagonal", "square"),
                        projection_crs = NULL,
                        id_column = NULL,
                        return_crs = 4326,
                        check_size = TRUE,
                        max_cells = 1000000) {
  if (!inherits(study_area, "sf")) stop("study_area must be sf", call. = FALSE)
  if (!is.numeric(cell_size) || length(cell_size) != 1L || cell_size <= 0) {
    stop("cell_size must be a single positive number", call. = FALSE)
  }
  type <- match.arg(type)

  if (is.null(projection_crs)) {
    projection_crs <- tryCatch(
      get_utm_crs(study_area),
      error = function(e) {
        stop(
          "Could not auto-detect a UTM projection_crs (", e$message,
          "). Pass projection_crs explicitly.",
          call. = FALSE
        )
      }
    )
    hex_msg("  Auto-selected projection_crs: ", projection_crs, "\n")
  }

  if (!is.null(id_column)) {
    if (!id_column %in% names(study_area)) {
      stop("id_column not found in study_area", call. = FALSE)
    }
    unique_areas <- unique(study_area[[id_column]])
    grids <- vector("list", length(unique_areas))
    names(grids) <- as.character(unique_areas)
    for (area_id in unique_areas) {
      hex_msg("Processing area: ", area_id, "\n")
      area_sf <- study_area[study_area[[id_column]] == area_id, ]
      grids[[as.character(area_id)]] <- create_single_grid(
        area_sf, cell_size, type, projection_crs, return_crs,
        area_id, check_size, max_cells
      )
    }
    return(grids)
  }

  create_single_grid(
    study_area, cell_size, type, projection_crs, return_crs,
    NULL, check_size, max_cells
  )
}

#' Find a hexagonal cell size that yields approximately a target number of cells
#'
#' Uses a binary search over `cell_size_min` / `cell_size_max` to find the
#' flat-to-flat distance that produces approximately `target_cells` hexagons.
#'
#' @param study_area sf object (polygons).
#' @param target_cells Desired number of hexagons (positive integer).
#' @param cell_size_min Minimum cell size to try. Default `NULL`, in which case
#'   the search range is derived from the bounding box of `study_area` so it
#'   works for either projected (metres) or geographic (degrees) input.
#' @param cell_size_max Maximum cell size to try. See `cell_size_min`.
#' @param tol Convergence tolerance (fraction of `target_cells`).
#' @param max_iter Maximum binary-search iterations.
#' @param projection_crs Optional CRS for grid construction (passed to
#'   [`create_grid()`]).
#' @return Cell size (flat-to-flat distance) closest to `target_cells`.
#' @export
find_hex_cell_size_for_target_cells <- function(study_area,
                                                target_cells,
                                                cell_size_min = NULL,
                                                cell_size_max = NULL,
                                                tol = 0.05,
                                                max_iter = 20,
                                                projection_crs = NULL) {
  if (!inherits(study_area, "sf")) stop("study_area must be sf", call. = FALSE)
  if (!is.numeric(target_cells) || length(target_cells) != 1L || target_cells <= 0) {
    stop("target_cells must be a single positive number", call. = FALSE)
  }

  if (is.null(cell_size_min) || is.null(cell_size_max)) {
    bb <- sf::st_bbox(study_area)
    diag <- sqrt((bb[3] - bb[1])^2 + (bb[4] - bb[2])^2)
    if (is.null(cell_size_min)) cell_size_min <- diag / 1e4
    if (is.null(cell_size_max)) cell_size_max <- diag / 2
  }

  lower <- cell_size_min
  upper <- cell_size_max
  best_size <- NA_real_
  best_diff <- Inf

  count_cells <- function(sz) {
    tryCatch({
      g <- create_grid(study_area, cell_size = sz, type = "hexagonal",
                       projection_crs = projection_crs, check_size = FALSE)
      nrow(g)
    }, error = function(e) NA_integer_)
  }

  hex_msg(sprintf(
    "[find_hex_cell_size_for_target_cells] Searching in [%.4f, %.4f]\n",
    lower, upper
  ))

  for (i in seq_len(max_iter)) {
    mid <- (lower + upper) / 2
    n_cells <- count_cells(mid)
    if (is.na(n_cells)) {
      lower <- mid
      next
    }
    diff_abs <- abs(n_cells - target_cells)
    if (diff_abs < best_diff) {
      best_diff <- diff_abs
      best_size <- mid
    }
    if (diff_abs / target_cells < tol) return(mid)
    if (n_cells > target_cells) {
      lower <- mid
    } else {
      upper <- mid
    }
  }

  if (is.na(best_size)) {
    warning(
      "Could not find a valid cell size in [", cell_size_min, ", ",
      cell_size_max, "]. Returning cell_size_max.",
      call. = FALSE
    )
    return(cell_size_max)
  }

  hex_msg(sprintf(
    "[find_hex_cell_size_for_target_cells] Closest: %.4f (diff = %d)\n",
    best_size, as.integer(best_diff)
  ))
  best_size
}

#' Build the grid for a single study area (internal)
#' @keywords internal
create_single_grid <- function(study_area, grid_size, type, projection_crs,
                               return_crs, area_id = NULL,
                               check_size = TRUE, max_cells = 1000000) {
  if (!inherits(study_area, "sf")) stop("study_area must be sf", call. = FALSE)
  if (nrow(study_area) == 0) stop("study_area is empty", call. = FALSE)
  if (grid_size <= 0) stop("grid_size must be positive", call. = FALSE)

  geom_types <- sf::st_geometry_type(study_area)
  if (!all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("study_area must contain (MULTI)POLYGON geometries", call. = FALSE)
  }

  study_area_proj <- tryCatch(
    sf::st_transform(study_area, projection_crs),
    error = function(e) stop("Failed to transform study_area: ", e$message, call. = FALSE)
  )
  bbox <- sf::st_bbox(study_area_proj)
  width <- bbox[3] - bbox[1]
  height <- bbox[4] - bbox[2]
  if (width <= 0 || height <= 0) stop("Invalid extent", call. = FALSE)

  hex_msg("  Creating ", type, " grid...\n")
  grid <- tryCatch(
    sf::st_sf(sf::st_make_grid(
      study_area_proj,
      cellsize = grid_size,
      square = (type == "square")
    )),
    error = function(e) {
      if (grepl("memory", e$message, ignore.case = TRUE)) {
        stop("Memory error: use a larger cell_size", call. = FALSE)
      }
      stop("Error making grid: ", e$message, call. = FALSE)
    }
  )

  if (check_size && nrow(grid) > 100000) {
    warning("Grid is very large (", nrow(grid), " cells)", call. = FALSE)
  }
  if (check_size && !is.null(max_cells) && nrow(grid) > max_cells) {
    stop(
      "Grid is too large (", nrow(grid), " cells). ",
      "Increase max_cells or use a larger cell_size.",
      call. = FALSE
    )
  }

  hex_msg("  Cropping grid...\n")
  grid_cropped <- tryCatch(
    sf::st_crop(grid, study_area_proj),
    error = function(e) stop("Error cropping: ", e$message, call. = FALSE)
  )
  rm(grid); gc(FALSE)

  hex_msg("  Clipping grid...\n")
  grid_clipped <- tryCatch(
    sf::st_intersection(grid_cropped, study_area_proj),
    error = function(e) stop("Error clipping: ", e$message, call. = FALSE)
  )
  rm(grid_cropped); gc(FALSE)

  hex_msg("  Transforming to output CRS...\n")
  grid_final <- tryCatch(
    sf::st_transform(grid_clipped, crs = return_crs),
    error = function(e) stop("Error transforming: ", e$message, call. = FALSE)
  )
  rm(grid_clipped); gc(FALSE)

  grid_final <- grid_final[
    sf::st_geometry_type(grid_final) %in% c("POLYGON", "MULTIPOLYGON"), ,
    drop = FALSE
  ]
  if (nrow(grid_final) == 0) stop("No valid grid cells created", call. = FALSE)

  if (!is.null(area_id)) {
    grid_final$area_id <- area_id
    grid_final$grid_id <- paste0(area_id, "_", seq_len(nrow(grid_final)))
  } else {
    grid_final$grid_id <- paste0("grid_", seq_len(nrow(grid_final)))
  }
  grid_final$grid_index <- seq_len(nrow(grid_final))

  hex_msg("  Grid created with ", nrow(grid_final), " cells\n")
  grid_final
}
