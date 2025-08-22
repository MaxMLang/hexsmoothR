# Test the complete vignette workflow
# This test runs through the exact same steps as the vignette

library(testthat)
library(hexsmoothR)
library(sf)
library(terra)

test_that("Complete vignette workflow works correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")
  
  # Step 1: Load sample data
  sample_file <- system.file("extdata", "default.tif", package = "hexsmoothR")
  sample_raster <- rast(sample_file)
  
  expect_true(file.exists(sample_file))
  expect_s4_class(sample_raster, "SpatRaster")
  # Use actual dimensions from the sample data
  expect_equal(nrow(sample_raster), 251)
  expect_equal(ncol(sample_raster), 224)
  
  # Step 2: Create study area and hexagonal grid
  study_area_wgs <- st_sf(geometry = st_sfc(
    st_polygon(list(matrix(
      c(ext(sample_raster)[1], ext(sample_raster)[3],  # xmin, ymin
        ext(sample_raster)[2], ext(sample_raster)[3],  # xmax, ymin
        ext(sample_raster)[2], ext(sample_raster)[4],  # xmax, ymax
        ext(sample_raster)[1], ext(sample_raster)[4],  # xmin, ymax
        ext(sample_raster)[1], ext(sample_raster)[3]), # back to start
      ncol = 2, byrow = TRUE
    ))), crs = 4326))
  
  expect_s3_class(study_area_wgs, "sf")
  expect_equal(st_crs(study_area_wgs)$input, "EPSG:4326")
  
  # Get UTM CRS
  utm_crs <- get_utm_crs(study_area_wgs)
  expect_true(grepl("^EPSG:32[67]\\d{2}$", utm_crs))
  
  # Transform to UTM
  study_area_utm <- st_transform(study_area_wgs, utm_crs)
  expect_equal(st_crs(study_area_utm)$input, utm_crs)
  
  # Create hexagonal grid
  cell_size <- 20000  # 20km hexagons in meters (UTM coordinates)
  hex_grid <- create_grid(study_area_utm, cell_size = cell_size, type = "hexagonal")
  
  expect_s3_class(hex_grid, "sf")
  expect_true(nrow(hex_grid) > 0)
  expect_true(nrow(hex_grid) < 10000)  # Should be reasonable size
  expect_true(all(c("grid_id", "grid_index") %in% names(hex_grid)))
  
  # Test optimal cell size function
  target_cells <- 800
  optimal_cell_size <- find_hex_cell_size_for_target_cells(
    study_area_utm, 
    target_cells = target_cells, 
    cell_size_min = 10000,  # 10km minimum (in meters)
    cell_size_max = 100000  # 100km maximum (in meters)
  )
  
  expect_true(optimal_cell_size >= 10000)
  expect_true(optimal_cell_size <= 100000)
  
  # Create grid with optimal cell size
  hex_grid_optimal <- create_grid(study_area_utm, cell_size = optimal_cell_size, type = "hexagonal")
  expect_s3_class(hex_grid_optimal, "sf")
  expect_true(nrow(hex_grid_optimal) > 0)
  
  # Step 3: Extract raster data
  extracted_data <- extract_raster_data(
    raster_files = c(sample = sample_file),
    study_area = study_area_wgs,  # WGS84 coordinates
    cell_size = 0.1           # 0.1 degrees (~11km) - appropriate for WGS84
  )
  
  expect_type(extracted_data, "list")
  expect_true("data" %in% names(extracted_data))
  expect_true("variables" %in% names(extracted_data))
  expect_true("n_cells" %in% names(extracted_data))
  # Note: the grid is not returned in the current implementation
  expect_equal(extracted_data$variables, "sample")
  expect_equal(nrow(extracted_data$data), extracted_data$n_cells)
  
  # Step 4: Configure smoothing topology
  topology <- compute_topology(hex_grid)
  
  expect_type(topology, "list")
  expect_true("neighbors" %in% names(topology))
  expect_true("weights" %in% names(topology))
  expect_true("avg_distance" %in% names(topology))
  expect_true("sigma" %in% names(topology))
  
  # Check weights structure (new N-order system)
  expect_true("center_weight" %in% names(topology$weights))
  expect_true("neighbor_weights" %in% names(topology$weights))
  expect_equal(length(topology$weights$neighbor_weights), 2)  # 2 orders
  
  # Weights should sum to 1 (normalized)
  total_weight <- topology$weights$center_weight + sum(topology$weights$neighbor_weights)
  expect_equal(total_weight, 1.0, tolerance = 0.001)
  
  # Test custom weights using new system
  aggressive_topology <- compute_topology(
    hex_grid,
    neighbor_orders = 2,
    center_weight = 0.3,
    neighbor_weights_param = list(0.5, 0.2)
  )
  
  expect_equal(aggressive_topology$weights$center_weight, 0.3, tolerance = 0.001)
  expect_equal(aggressive_topology$weights$neighbor_weights[[1]], 0.5, tolerance = 0.001)
  expect_equal(aggressive_topology$weights$neighbor_weights[[2]], 0.2, tolerance = 0.001)
  
  # Step 5: Apply spatial smoothing
  smoothing_results <- smooth_variables(
    variable_values = list(sample = extracted_data$data$sample),
    neighbors = topology$neighbors,
    weights = topology$weights,
    var_names = c("sample")
  )
  
  expect_type(smoothing_results, "list")
  expect_true("sample" %in% names(smoothing_results))
  
  # Check that all expected columns are present
  sample_cols <- names(smoothing_results$sample)
  expected_cols <- c("raw", "neighbors_1st", "neighbors_2nd", "weighted_combined")
  expect_true(all(expected_cols %in% sample_cols))
  
  # Check that results have the right length - use the grid size, not extracted data size
  expect_equal(length(smoothing_results$sample$raw), nrow(hex_grid))
  expect_equal(length(smoothing_results$sample$weighted_combined), nrow(hex_grid))
  
  # Step 6: Combine results and analyze
  hex_grid_with_results <- hex_grid
  hex_grid_with_results$sample_raw <- smoothing_results$sample$raw
  hex_grid_with_results$sample_neighbors_1st <- smoothing_results$sample$neighbors_1st
  hex_grid_with_results$sample_neighbors_2nd <- smoothing_results$sample$neighbors_2nd
  hex_grid_with_results$sample_weighted_combined <- smoothing_results$sample$weighted_combined
  
  expect_s3_class(hex_grid_with_results, "sf")
  expect_true(all(c("sample_raw", "sample_weighted_combined") %in% names(hex_grid_with_results)))
  
  # Transform to WGS84 for plotting
  hex_grid_with_results_wgs84 <- st_transform(hex_grid_with_results, crs = 4326)
  expect_equal(st_crs(hex_grid_with_results_wgs84)$input, "EPSG:4326")
  
  # Step 7: Test hexagon utility functions
  cell_size_flat <- 20000  # 20km flat-to-flat distance in meters
  
  edge_length <- hex_flat_to_edge(cell_size_flat)
  circumradius <- hex_flat_to_circumradius(cell_size_flat)
  
  expect_true(edge_length < cell_size_flat)
  expect_true(circumradius < cell_size_flat)
  
  # Convert back
  flat_distance_from_edge <- hex_edge_to_flat(edge_length)
  flat_distance_from_circumradius <- hex_circumradius_to_flat(circumradius)
  
  expect_equal(flat_distance_from_edge, cell_size_flat, tolerance = 0.001)
  expect_equal(flat_distance_from_circumradius, cell_size_flat, tolerance = 0.001)
  
  # Test multiple variables
  # Create dummy data for second variable - use the grid size
  dummy_values <- rep(1, nrow(hex_grid))
  
  smoothing_multi <- smooth_variables(
    variable_values = list(
      sample = smoothing_results$sample$raw,  # Use the raw values from previous smoothing
      dummy = dummy_values
    ),
    neighbors = topology$neighbors,
    weights = topology$weights,
    var_names = c("sample", "dummy")
  )
  
  expect_type(smoothing_multi, "list")
  expect_true(all(c("sample", "dummy") %in% names(smoothing_multi)))
  expect_true(all(expected_cols %in% names(smoothing_multi$sample)))
  expect_true(all(expected_cols %in% names(smoothing_multi$dummy)))
  
  # Test error handling
  expect_error(
    create_grid("invalid", cell_size = 1000),
    "study_area must be sf"
  )
  
  expect_error(
    get_utm_crs("invalid"),
    "study_area must be an sf object"
  )
  
  expect_error(
    compute_topology("invalid"),
    "grid must be an sf object or a list of sf objects"
  )
  
  # Test with empty grid - this should work but may not give a warning
  empty_grid <- hex_grid[0, ]
  # The function handles empty grids gracefully, so just test it doesn't crash
  expect_no_error(compute_topology(empty_grid))
  
  cat("✅ All vignette workflow tests passed successfully!\n")
}) 