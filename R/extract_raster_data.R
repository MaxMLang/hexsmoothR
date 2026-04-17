#' Extract raster data into a hexagonal grid
#'
#' Primary function for extracting raster values into hexagonal grid cells. CRS
#' transformations are handled automatically: each raster is reprojected (or the
#' grid is reprojected to the raster's CRS) so that the underlying call to
#' [`exactextractr::exact_extract()`] always sees matching CRSs.
#'
#' @section Input flexibility:
#' `raster_files` may be either a named character vector of file paths or a
#' named list of `terra::SpatRaster` objects. All inputs may be in different
#' CRSs from one another and from the grid - the function handles cropping and
#' transformation per raster.
#'
#' @section CRS handling:
#' - If `study_area` is supplied, the grid is created in the study area's CRS.
#' - Otherwise the grid is created in the first raster's CRS.
#' - For each raster, the grid is transformed to the raster's CRS before
#'   extraction (so cell-size units are honoured exactly once, in the grid CRS).
#'
#' @param raster_files Named character vector of file paths OR named list of
#'   `terra::SpatRaster` objects.
#' @param study_area Optional sf polygon used for cropping each raster and to
#'   define the grid CRS.
#' @param cell_size Hex cell size, in the units of the grid CRS. Required when
#'   `hex_grid` is not supplied.
#' @param hex_grid Optional sf hexagonal grid to use instead of creating one.
#' @param sample_fraction Fraction of grid cells to keep (default 1).
#' @param random_seed Seed for reproducible sampling.
#' @param fun Aggregation function passed to `exactextractr::exact_extract()`
#'   (default `"mean"`).
#' @return List with components
#'   \describe{
#'     \item{`data`}{Data frame with `cell_id`, `x`, `y` and one column per raster.}
#'     \item{`hex_grid`}{The sf grid that was used (sampled, if applicable).}
#'     \item{`cell_size`}{The cell size used.}
#'     \item{`extent`}{Extent of the first raster (after cropping).}
#'     \item{`variables`}{Names of the rasters.}
#'     \item{`n_cells`}{Number of cells in `data`.}
#'   }
#'
#' @importFrom terra rast res vect crop ncell ext crs
#' @importFrom sf st_transform st_crs st_make_grid st_geometry st_centroid
#'   st_coordinates st_is_longlat st_sf
#' @importFrom exactextractr exact_extract
#' @export
extract_raster_data <- function(raster_files,
                                study_area = NULL,
                                cell_size = NULL,
                                hex_grid = NULL,
                                sample_fraction = 1.0,
                                random_seed = 42,
                                fun = "mean") {
  if (is.null(names(raster_files)) || any(!nzchar(names(raster_files)))) {
    stop("raster_files must be fully named", call. = FALSE)
  }
  if (sample_fraction <= 0 || sample_fraction > 1) {
    stop("sample_fraction must be in (0, 1]", call. = FALSE)
  }
  set.seed(random_seed)

  is_file_paths  <- all(vapply(raster_files, is.character, logical(1)))
  is_rast_objects <- all(vapply(raster_files, inherits, logical(1), "SpatRaster"))
  if (!is_file_paths && !is_rast_objects) {
    stop(
      "raster_files must be either all character file paths or all ",
      "terra::SpatRaster objects",
      call. = FALSE
    )
  }

  load_raster <- function(idx) {
    if (is_file_paths) terra::rast(raster_files[[idx]]) else raster_files[[idx]]
  }

  hex_msg(
    "Extracting raster data from ", length(raster_files), " ",
    if (is_file_paths) "files" else "SpatRaster objects", "...\n"
  )

  first_raster <- load_raster(1)
  raster_crs <- terra::crs(first_raster)
  hex_msg(
    "First raster: ", names(raster_files)[1],
    " (", nrow(first_raster), " x ", ncol(first_raster), ")\n"
  )

  if (!is.null(study_area)) {
    if (!inherits(study_area, "sf")) stop("study_area must be sf", call. = FALSE)
    grid_crs <- sf::st_crs(study_area)
  } else {
    grid_crs <- sf::st_crs(raster_crs)
  }
  is_grid_geographic <- isTRUE(sf::st_is_longlat(grid_crs))

  if (!is.null(cell_size) && is_grid_geographic && cell_size > 1) {
    warning(
      "cell_size = ", cell_size, " is large for a geographic CRS (degrees). ",
      "Consider a projected CRS or a smaller cell_size.",
      call. = FALSE
    )
  }

  if (!is.null(study_area)) {
    study_area_for_first <- sf::st_transform(study_area, sf::st_crs(raster_crs))
    first_raster <- terra::crop(first_raster, terra::vect(study_area_for_first))
    hex_msg("Cropped first raster to study area\n")
  }

  if (!is.null(hex_grid)) {
    if (!inherits(hex_grid, "sf")) stop("hex_grid must be sf", call. = FALSE)
    hex_grid_sf <- hex_grid
  } else {
    if (is.null(cell_size)) {
      cell_size <- max(terra::res(first_raster)) * 10
      hex_msg("Using automatic cell_size: ", cell_size, "\n")
    }
    grid_template <- if (!is.null(study_area)) study_area else sf::st_as_sf(
      sf::st_as_sfc(sf::st_bbox(first_raster), crs = sf::st_crs(raster_crs))
    )
    hex_msg("Creating hexagonal grid in grid CRS...\n")
    hex_geom <- sf::st_make_grid(
      grid_template,
      cellsize = cell_size,
      square = FALSE,
      what = "polygons"
    )
    hex_grid_sf <- sf::st_sf(geometry = hex_geom, crs = sf::st_crs(grid_template))
  }

  if (!"cell_id" %in% names(hex_grid_sf)) {
    hex_grid_sf$cell_id <- seq_len(nrow(hex_grid_sf))
  }

  if (sample_fraction < 1) {
    n_total <- nrow(hex_grid_sf)
    n_samples <- round(n_total * sample_fraction)
    sample_indices <- sample(seq_len(n_total), n_samples)
    hex_grid_sf <- hex_grid_sf[sample_indices, , drop = FALSE]
    hex_msg("Sampled ", n_samples, " cells from ", n_total, " total cells\n")
  }

  centroids <- suppressWarnings(sf::st_centroid(hex_grid_sf))
  coords <- sf::st_coordinates(centroids)

  result_data <- data.frame(
    cell_id = hex_grid_sf$cell_id,
    x = coords[, 1],
    y = coords[, 2],
    stringsAsFactors = FALSE
  )

  for (var_name in names(raster_files)) {
    hex_msg("Processing ", var_name, "...\n")

    r <- tryCatch(
      load_raster(var_name),
      error = function(e) {
        hex_msg("  Error loading ", var_name, ": ", e$message, "\n")
        NULL
      }
    )
    if (is.null(r) || !inherits(r, "SpatRaster")) {
      result_data[[var_name]] <- NA_real_
      next
    }

    if (!is.null(study_area)) {
      sa_in_r_crs <- sf::st_transform(study_area, sf::st_crs(terra::crs(r)))
      r <- terra::crop(r, terra::vect(sa_in_r_crs))
    }

    grid_in_r_crs <- if (sf::st_crs(hex_grid_sf) == sf::st_crs(terra::crs(r))) {
      hex_grid_sf
    } else {
      sf::st_transform(hex_grid_sf, sf::st_crs(terra::crs(r)))
    }

    extracted <- tryCatch(
      exactextractr::exact_extract(r, grid_in_r_crs, fun, progress = FALSE),
      error = function(e) {
        hex_msg("  Error extracting ", var_name, ": ", e$message, "\n")
        rep(NA_real_, nrow(hex_grid_sf))
      }
    )
    result_data[[var_name]] <- extracted
    hex_msg(
      "  Extracted ", sum(!is.na(extracted)), " valid values out of ",
      length(extracted), "\n"
    )
  }

  result <- list(
    data = result_data,
    hex_grid = hex_grid_sf,
    cell_size = cell_size,
    extent = terra::ext(first_raster),
    variables = names(raster_files),
    n_cells = nrow(result_data)
  )

  hex_msg(
    "Extraction done. Cells: ", nrow(result_data),
    "; variables: ", paste(names(raster_files), collapse = ", "), "\n"
  )
  result
}
