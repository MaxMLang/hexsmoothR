#' Extract raster data for a hex grid using exactextractr
#'
#' This is the PRIMARY function for extracting raster data to hexagonal grids.
#' It automatically handles CRS transformations and creates consistent hexagonal grids.
#'
#' **CRS Handling:**
#' - If `study_area` is provided, the function creates a grid in the study area's CRS
#' - If no `study_area` is provided, the function creates a grid in the raster's CRS
#' - The function automatically crops rasters to the study area (if provided)
#' - All coordinate transformations are handled transparently
#'
#' **Grid Creation:**
#' - Creates hexagonal grids using `sf::st_make_grid` for robustness
#' - Automatically adjusts cell sizes for geographic coordinates (degrees)
#' - Can use a pre-existing grid via the `hex_grid` parameter
#'
#' **Input Flexibility:**
#' - Accepts both file paths (character strings) and `terra::rast` objects
#' - When passing `terra::rast` objects, you can inspect CRS, dimensions, and data before extraction
#' - Useful for working with in-memory rasters or pre-processed data
#'
#' @param raster_files named vector of raster file paths OR named list of terra::rast objects
#' @param study_area sf polygon (optional, for cropping and grid CRS reference)
#' @param cell_size hex cell size (in units of the target CRS - see details)
#' @param hex_grid optional sf hexagonal grid to use instead of creating a new one
#' @param sample_fraction fraction of raster cells to use (default 1)
#' @param random_seed for reproducibility
#' @return list: data frame with cell_id, x, y, and values for each raster
#' 
#' @details
#' **CRS Handling Details:**
#' 
#' The function intelligently handles coordinate reference systems:
#' 
#' 1. **Study Area Provided**: Grid is created in the study area's CRS
#'    - Example: If study_area is in UTM (EPSG:32630), grid uses UTM coordinates
#'    - Cell sizes should be in meters (e.g., 20000 for 20km hexagons)
#' 
#' 2. **No Study Area**: Grid is created in the raster's CRS
#'    - Example: If raster is in WGS84 (EPSG:4326), grid uses geographic coordinates
#'    - Cell sizes should be in degrees (e.g., 0.1 for ~11km hexagons)
#' 
#' 3. **Automatic Adjustments**: 
#'    - Large cell sizes (>1) in geographic coordinates trigger warnings
#'    - Cell sizes are automatically adjusted for geographic coordinates if needed
#' 
#' **Input Types:**
#' 
#' The function accepts two types of input for `raster_files`:
#' 
#' 1. **File Paths**: Named character vector of file paths
#'    ```r
#'    raster_files = c(ndvi = "path/to/ndvi.tif", elevation = "path/to/elevation.tif")
#'    ```
#' 
#' 2. **Terra Raster Objects**: Named list of terra::rast objects
#'    ```r
#'    # Load rasters first to inspect them
#'    ndvi_rast <- rast("path/to/ndvi.tif")
#'    elev_rast <- rast("path/to/elevation.tif")
#'    
#'    # Check CRS, dimensions, etc.
#'    print(crs(ndvi_rast))
#'    print(dim(ndvi_rast))
#'    
#'    # Pass to function
#'    raster_files = list(ndvi = ndvi_rast, elevation = elev_rast)
#'    ```
#' 
#' **Workflow Recommendation:**
#' 
#' For consistent results, always provide a `study_area` in your desired CRS:
#' ```r
#' # Create study area in UTM for accurate measurements
#' study_area_utm <- st_transform(study_area_wgs, utm_crs)
#' 
#' # Extract with 20km hexagons (in UTM coordinates)
#' extracted <- extract_raster_data(
#'   raster_files = c(var1 = "path/to/raster.tif"),
#'   study_area = study_area_utm,
#'   cell_size = 20000  # 20km in meters
#' )
#' 
#' # Then apply smoothing
#' smoothed <- smooth_variables(...)
#' ```
#' 
#' @importFrom terra rast res vect crop ncell ext crs
#' @importFrom sf st_transform st_crs st_make_grid st_geometry
#' @importFrom exactextractr exact_extract
#' @export
extract_raster_data <- function(raster_files, 
                              study_area = NULL,
                              cell_size = NULL,
                              hex_grid = NULL,
                              sample_fraction = 1.0,
                              random_seed = 42) {
  if (is.null(names(raster_files))) stop("raster_files must be named")
  if (sample_fraction <= 0 || sample_fraction > 1) stop("sample_fraction must be 0-1")
  set.seed(random_seed)
  
  # Determine input type and process accordingly
  is_file_paths <- all(sapply(raster_files, is.character))
  is_rast_objects <- all(sapply(raster_files, function(x) inherits(x, "SpatRaster")))
  
  if (!is_file_paths && !is_rast_objects) {
    stop("raster_files must be either all file paths (character) or all terra::rast objects")
  }
  
  if (is_file_paths) {
    cat("Extracting raster data from", length(raster_files), "files...\n")
  } else {
    cat("Extracting raster data from", length(raster_files), "terra::rast objects...\n")
  }
  
  # Load or use first raster to get reference CRS and properties
  if (is_file_paths) {
    first_raster <- rast(raster_files[1])
    cat("First raster loaded:", names(raster_files)[1], "\n")
  } else {
    first_raster <- raster_files[[1]]
    cat("First raster object:", names(raster_files)[1], "\n")
  }
  
  cat("Dimensions:", nrow(first_raster), "x", ncol(first_raster), "\n")
  cat("Resolution:", res(first_raster), "\n")
  cat("CRS:", crs(first_raster), "\n")
  
  # Determine the target coordinate system for grid creation
  raster_crs <- crs(first_raster)
  
  # CRS handling logic - be more explicit about what we're doing
  if (!is.null(study_area)) {
    study_area_crs <- st_crs(study_area)$wkt
    
    # Check if study area and raster have different CRS
    if (study_area_crs != raster_crs) {
      cat("Study area and raster have different CRS.\n")
      cat("Study area CRS:", substr(study_area_crs, 1, 80), "...\n")
      cat("Raster CRS:", substr(raster_crs, 1, 80), "...\n")
      
      # Determine if study area is projected (UTM) and raster is geographic
      is_study_projected <- !grepl("GEOGCRS|longlat|latlong", study_area_crs)
      is_raster_geographic <- grepl("GEOGCRS|longlat|latlong", raster_crs)
      
      if (is_study_projected && is_raster_geographic) {
        cat("Study area is projected (UTM) but raster is geographic (WGS84).\n")
        cat("This is the recommended setup for accurate measurements.\n")
        cat("Grid will be created in study area CRS (UTM) with cell_size in meters.\n")
        target_crs <- study_area_crs
        grid_crs <- study_area_crs
      } else if (is_study_projected && !is_raster_geographic) {
        cat("Both study area and raster are projected. Using study area CRS.\n")
        target_crs <- study_area_crs
        grid_crs <- study_area_crs
      } else {
        cat("Using raster CRS for grid creation.\n")
        target_crs <- raster_crs
        grid_crs <- raster_crs
      }
    } else {
      cat("Study area and raster have same CRS. Using this CRS.\n")
      target_crs <- raster_crs
      grid_crs <- raster_crs
    }
  } else {
    cat("No study area provided. Using raster CRS for grid creation.\n")
    target_crs <- raster_crs
    grid_crs <- raster_crs
  }
  
  # Cell size validation and adjustment
  if (!is.null(cell_size)) {
    # Check if we're working with geographic coordinates
    is_grid_geographic <- grepl("GEOGCRS|longlat|latlong", grid_crs)
    
    if (is_grid_geographic) {
      # Grid will be in geographic coordinates (degrees)
      if (cell_size > 1) {
        cat("Warning: Large cell_size (", cell_size, ") detected for geographic coordinates.\n")
        cat("This may cause issues. Consider using a smaller cell_size (< 1 degree) or projected coordinates.\n")
        
        # Suggest a reasonable cell size for degrees
        suggested_cell_size <- 0.1  # 0.1 degrees is roughly 11km at the equator
        cat("Suggested cell_size for degrees:", suggested_cell_size, "\n")
        
        # Only adjust if the cell_size is unreasonably large
        if (cell_size > 10) {
          cell_size <- suggested_cell_size
          cat("Adjusted cell_size to:", cell_size, "degrees\n")
        }
      }
    } else {
      # Grid will be in projected coordinates (meters)
      cat("Grid will be created in projected coordinates (UTM).\n")
      cat("Cell size should be in meters (e.g., 20000 for 20km hexagons).\n")
    }
  }
  
  # Handle study area and raster cropping
  if (!is.null(study_area)) {
    study_area_crs <- st_crs(study_area)$wkt
    
    if (raster_crs != study_area_crs) {
      # Transform study area to raster CRS for cropping only
      cat("Transforming study area to raster CRS for cropping...\n")
      study_area_for_crop <- st_transform(study_area, raster_crs)
      study_area_vect <- vect(study_area_for_crop)
      first_raster <- crop(first_raster, study_area_vect)
      cat("Cropped to study area\n")
    } else {
      # Same CRS, no transformation needed
      study_area_vect <- vect(study_area)
      first_raster <- crop(first_raster, study_area_vect)
      cat("Cropped to study area\n")
    }
  }
  
  # Use provided grid or create new one
  if (!is.null(hex_grid)) {
    cat("Using provided hexagonal grid...\n")
    hex_grid_sf <- hex_grid
    
    # Transform grid to raster CRS if needed
    if (st_crs(hex_grid_sf)$wkt != raster_crs) {
      cat("Transforming grid to raster CRS...\n")
      hex_grid_sf <- st_transform(hex_grid_sf, raster_crs)
    }
    
    # Add cell IDs if missing
    if (!"cell_id" %in% names(hex_grid_sf)) {
      hex_grid_sf$cell_id <- 1:nrow(hex_grid_sf)
    }
    
  } else {
    # Create new grid
    # Determine cell size (if not provided)
    if (is.null(cell_size)) {
      cell_size <- max(res(first_raster)) * 10
      cat("Using automatic cell size:", cell_size, "\n")
    }
    
    # Create hexagonal grid using sf::st_make_grid (more robust than custom hex generation)
    cat("Creating hexagonal grid...\n")
    
    # Always create grid in the target CRS
    if (grid_crs != raster_crs) {
      # Create grid in the target CRS using the study area
      cat("Creating grid in target CRS (", substr(grid_crs, 1, 50), "...)\n")
      hex_grid <- st_make_grid(
        study_area, 
        cellsize = cell_size, 
        square = FALSE,  # This creates hexagons
        what = "polygons"
      )
      hex_grid_sf <- sf::st_sf(geometry = hex_grid, crs = grid_crs)
    } else {
      # Create grid in the raster's CRS
      cat("Creating grid in raster CRS\n")
      hex_grid <- st_make_grid(
        first_raster, 
        cellsize = cell_size, 
        square = FALSE,  # This creates hexagons
        what = "polygons"
      )
      hex_grid_sf <- sf::st_sf(geometry = hex_grid, crs = raster_crs)
    }
    
    # Add cell IDs
    hex_grid_sf$cell_id <- 1:nrow(hex_grid_sf)
  }
  
  # Sample grid if needed
  if (sample_fraction < 1) {
    n_samples <- round(nrow(hex_grid_sf) * sample_fraction)
    sample_indices <- sample(1:nrow(hex_grid_sf), n_samples)
    hex_grid_sf <- hex_grid_sf[sample_indices, ]
    cat("Sampled", n_samples, "cells from", nrow(hex_grid_sf), "total cells\n")
  }
  
  # Get centroids for coordinate output
  centroids <- st_centroid(hex_grid_sf)
  coords <- st_coordinates(centroids)
  
  # Initialize result data frame
  result_data <- data.frame(
    cell_id = hex_grid_sf$cell_id,
    x = coords[, 1], 
    y = coords[, 2], 
    stringsAsFactors = FALSE
  )
  
  # Extract values from each raster using exactextractr
  for (var_name in names(raster_files)) {
    cat("Processing", var_name, "...\n")
    
    # Load or use raster based on input type
    if (is_file_paths) {
      r <- tryCatch(
        rast(raster_files[var_name]), 
        error = function(e) { 
          cat("  Error loading", var_name, ":", e$message, "\n"); 
          return(NULL) 
        }
      )
    } else {
      r <- raster_files[[var_name]]
      if (!inherits(r, "SpatRaster")) {
        cat("  Error: ", var_name, "is not a terra::rast object\n")
        result_data[[var_name]] <- rep(NA, nrow(hex_grid_sf))
        next
      }
    }
    
    if (is.null(r)) { 
      result_data[[var_name]] <- NA
      next 
    }
    
    # Crop raster if study area provided
    if (!is.null(study_area)) {
      r <- crop(r, study_area_vect)
    }
    
    # Extract values using exactextractr (much more robust than terra::extract)
    tryCatch({
      values <- exact_extract(r, hex_grid_sf, 'mean')
      result_data[[var_name]] <- values
      cat("  Extracted", sum(!is.na(values)), "valid values out of", length(values), "total\n")
    }, error = function(e) {
      cat("  Error extracting values for", var_name, ":", e$message, "\n")
      result_data[[var_name]] <<- rep(NA, nrow(hex_grid_sf))
    })
  }
  
  result <- list(
    data = result_data, 
    cell_size = cell_size, 
    extent = ext(first_raster), 
    variables = names(raster_files), 
    n_cells = nrow(result_data)
  )
  
  cat("Extraction done!\nTotal cells:", nrow(result_data), "\nVariables:", paste(names(raster_files), collapse = ", "), "\n\n")
  return(result)
}

 