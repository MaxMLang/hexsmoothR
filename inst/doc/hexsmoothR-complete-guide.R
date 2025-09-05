## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10,
  fig.height = 7,
  warning = FALSE,
  message = FALSE
)

## ----load-libraries, eval = TRUE----------------------------------------------
library(hexsmoothR)
library(sf)
library(terra)
library(exactextractr)

# Define area of interest for coastal analysis
bbox <- st_bbox(
  c(xmin = -122.50, ymin = 47.62, xmax = -122.54, ymax = 47.63),
  crs = st_crs(4326)
)

aoi_sf <- st_as_sf(st_as_sfc(bbox))
cat("AOI: Puget Sound coastal region\n")

# Load Sentinel-2 L2A NDVI data or fallback to simulated data
sentinel_file <- system.file("extdata", "sentinel2_ndvi.tif", package = "hexsmoothR")
if (file.exists(sentinel_file)) {
  true_coastal_vegetation <- rast(sentinel_file) 
  cat("Sentinel-2 L2A NDVI data loaded\n")
} else {
  cat("Using simulated coastal vegetation data\n")
  
  # Create a simple simulated elevation raster for demonstration
  set.seed(42)
  true_elevation <- rast(nrows = 100, ncols = 100, 
                         xmin = -122.54, xmax = -122.50, 
                         ymin = 47.62, ymax = 47.63, 
                         crs = "EPSG:4326")
  
  # Simulate elevation with coastal gradient
  x_coords <- seq(-122.54, -122.50, length.out = 100)
  y_coords <- seq(47.62, 47.63, length.out = 100)
  coords <- expand.grid(x = x_coords, y = y_coords)
  
  # Create elevation pattern: higher inland, lower near coast
  elevation_values <- 50 + 100 * (coords$x - min(coords$x)) / (max(coords$x) - min(coords$x)) +
                     30 * (coords$y - min(coords$y)) / (max(coords$y) - min(coords$y)) +
                     rnorm(nrow(coords), 0, 10)
  
  values(true_elevation) <- elevation_values
  
  # Simulate coastal vegetation pattern based on elevation
  vegetation_from_elev <- 1 - (true_elevation / 200)
  distance_effect <- (true_elevation - min(values(true_elevation), na.rm = TRUE)) / 
                    (max(values(true_elevation), na.rm = TRUE) - min(values(true_elevation), na.rm = TRUE))
  vegetation_from_distance <- 1 - distance_effect * 0.6
  
  set.seed(42)
  coastal_patches <- rast(true_elevation)
  values(coastal_patches) <- runif(ncell(true_elevation), 0, 0.3)
  
  true_coastal_vegetation <- vegetation_from_elev + vegetation_from_distance + coastal_patches
  true_coastal_vegetation[true_coastal_vegetation < 0.1] <- 0.1
  true_coastal_vegetation[true_coastal_vegetation > 0.8] <- 0.8
}

# Use original data for hexagonal smoothing demonstration
observed_coastal_vegetation <- true_coastal_vegetation

# Data summary
if (file.exists(sentinel_file)) {
  cat("Sentinel-2 L2A NDVI data loaded\n")
} else {
  cat("Simulated coastal vegetation data\n")
}
cat("Dimensions:", nrow(observed_coastal_vegetation), "x", ncol(observed_coastal_vegetation), "pixels\n")
cat("Resolution:", round(xres(observed_coastal_vegetation) * 111320, 1), "m\n")
cat("CRS:", crs(observed_coastal_vegetation, proj = TRUE), "\n")
missing_pct <- round(global(is.na(observed_coastal_vegetation), "mean")$mean * 100, 1)
cat("Missing values:", missing_pct, "%\n")

# Create color palette for NDVI visualization
plot_colors <- colorRampPalette(c("darkred", "red", "orange", "yellow", "green", "darkgreen"))(100)
cat("Color palette: 100 colors (red to green)\n")

# Plot 1: True underlying coastal vegetation pattern (NDVI-like)
if (file.exists(sentinel_file)) {
  plot_title_1 <- "Real Sentinel-2 L2A Coastal Vegetation Pattern (NDVI)\nThe signal to recover"
} else {
  plot_title_1 <- "True Coastal Vegetation Pattern (NDVI-like)\nThe signal to recover"
}

plot(true_coastal_vegetation$sentinel2_ndvi_1, 
     main = plot_title_1,
     col = plot_colors,
     axes = FALSE)

# Plot 2: Observed coastal vegetation data
if (file.exists(sentinel_file)) {
  plot_title_2 <- "Observed Coastal Vegetation Data (NDVI)\nOriginal Sentinel-2 data"
} else {
  plot_title_2 <- "Observed Coastal Vegetation Data (NDVI-like)\nSimulated data"
}

plot(observed_coastal_vegetation$sentinel2_ndvi_3,
     main = plot_title_2,
     col = plot_colors,
     axes = FALSE)

cat("Individual plots created successfully!\n")

demo_raster <- observed_coastal_vegetation

## ----create-grid, eval = TRUE-------------------------------------------------
# Create study area polygon from raster extent
study_area_wgs <- st_sf(geometry = st_sfc(
  st_polygon(list(matrix(
    c(ext(demo_raster)[1], ext(demo_raster)[3],  # xmin, ymin
      ext(demo_raster)[2], ext(demo_raster)[3],  # xmax, ymin
      ext(demo_raster)[2], ext(demo_raster)[4],  # xmax, ymax
      ext(demo_raster)[1], ext(demo_raster)[4],  # xmin, ymax
      ext(demo_raster)[1], ext(demo_raster)[3]), # back to start
    ncol = 2, byrow = TRUE
  ))), crs = 4326))

cat("Study area CRS: WGS84\n")

# Transform to UTM for accurate measurements
utm_crs <- get_utm_crs(study_area_wgs)
study_area_utm <- st_transform(study_area_wgs, utm_crs)

cat("UTM CRS:", st_crs(study_area_utm)$input, "\n")

# Calculate cell size for target hexagon count
# Convert degrees to meters (111320 factor, varies by latitude)
raster_area_km2 <- (nrow(demo_raster) * ncol(demo_raster) * (xres(demo_raster) * 111320/1000) * (yres(demo_raster) * 111320/1000))

# Target hexagon count for spatial resolution
target_cells <- 5000

# Calculate cell size using hexagon area ratio (2.598)
cell_size_meters <- sqrt((raster_area_km2 * 1e6) / (target_cells * 2.598))

cat("Raster area:", round(raster_area_km2, 1), "km²\n")
cat("Target hexagons:", target_cells, "\n")
cat("Cell size:", round(cell_size_meters, 0), "m (", round(cell_size_meters/1000, 1), "km)\n")

hex_grid <- create_grid(study_area_utm, cell_size = cell_size_meters, type = "hexagonal")

cat("Grid created:", nrow(hex_grid), "cells\n")
cat("Average cell area:", round(mean(st_area(hex_grid)) / 1e6, 1), "km²\n")

# Transform to WGS84 for raster compatibility
hex_grid_wgs84 <- st_transform(hex_grid, crs = 4326)

# Visualize grid overlay
par(mfrow = c(1, 2))

plot(demo_raster$sentinel2_ndvi_1, main = "Data with Hexagonal Grid")
plot(st_geometry(hex_grid_wgs84), add = TRUE, border = "white", lwd = 0.1)

plot(st_geometry(hex_grid_wgs84), main = paste("Hexagonal Grid -", nrow(hex_grid), "cells"), 
     border = "blue", lwd = 0.05)

## ----optimal-cell-size, eval = TRUE-------------------------------------------
# Package suggests optimal cell size for target hexagon count
target_cells <- 5000
optimal_cell_size <- find_hex_cell_size_for_target_cells(
  study_area_utm, 
  target_cells = target_cells, 
  cell_size_min = 20,  
  cell_size_max = 40
)

cat("Optimal cell size:", round(optimal_cell_size, 0), "m\n")

## ----import-grid, eval = FALSE------------------------------------------------
# # Load existing hexagonal grid and ensure correct CRS
# existing_grid <- st_read("path/to/your/hex_grid.shp")
# 
# if (st_crs(existing_grid)$input != st_crs(study_area_utm)$input) {
#   existing_grid <- st_transform(existing_grid, crs = st_crs(study_area_utm))
# }
# 
# hex_grid <- existing_grid

## ----extract-raster, eval = TRUE----------------------------------------------
# Extract raster values into hexagon cells
# Function handles CRS transformations and pixel aggregation

cat("Raster data extraction\n")
cat("Using existing grid (cell_size = NULL)\n")

# Extract raster values using existing hexagonal grid
extracted_data <- extract_raster_data(
  raster_files = list(coastal_vegetation = demo_raster),
  hex_grid = hex_grid_wgs84,
  cell_size = NULL
)

cat("Extraction complete:", extracted_data$n_cells, "cells\n")
cat("Variables:", paste(extracted_data$variables, collapse = ", "), "\n")
cat("Valid values:", sum(!is.na(extracted_data$data$coastal_vegetation)), "/", nrow(extracted_data$data), "\n")

# Data statistics comparison
original_mean <- mean(values(demo_raster), na.rm = TRUE)
original_sd <- sd(values(demo_raster), na.rm = TRUE)

# Handle data structure
if (is.list(extracted_data$data$coastal_vegetation)) {
  extracted_values <- extracted_data$data$coastal_vegetation[[1]]
} else {
  extracted_values <- extracted_data$data$coastal_vegetation
}

extracted_mean <- mean(extracted_values, na.rm = TRUE)
extracted_sd <- sd(extracted_values, na.rm = TRUE)

cat("Data comparison:\n")
cat("Original raster - Mean:", round(original_mean, 3), "SD:", round(original_sd, 3), "\n")
cat("Extracted hexagons - Mean:", round(extracted_mean, 3), "SD:", round(extracted_sd, 3), "\n")

# Visualize extraction results
plot(demo_raster$sentinel2_ndvi_1, 
     main = paste("Original Raster (", ncell(demo_raster), "pixels)"),
     col = plot_colors)

# Plot extracted hexagon values
hex_for_plot <- hex_grid_wgs84
hex_for_plot$value <- extracted_values

# Use base R plotting for consistency
plot(st_geometry(hex_for_plot), 
     col = plot_colors[cut(hex_for_plot$value, 100)],
     main = paste("Extracted Hexagon Coastal Vegetation Values (NDVI-like)\n(", extracted_data$n_cells, "hexagons)"),
     border = "white", lwd = 0.01)

cat("\nExtraction process:\n")
cat("- Left: Original coastal vegetation raster with", ncell(demo_raster), "pixels and noise\n")  
cat("- Right: Extracted coastal vegetation values in", extracted_data$n_cells, "hexagons\n")
cat("- Each hexagon contains the average coastal vegetation of all pixels within its boundary\n")
cat("- This aggregation already reduces some measurement noise!\n")
cat("- **PERFECT ALIGNMENT**: extract_raster_data preserves grid orientation like exact_extract\n")

# The extracted data contains:
# - data: dataframe with extracted values
# - variables: names of the extracted variables
# - n_cells: number of grid cells

## ----multiple-rasters, eval = FALSE-------------------------------------------
# # Load multiple rasters for inspection
# ndvi_raster <- rast("path/to/ndvi.tif")
# elevation_raster <- rast("path/to/elevation.tif")
# precipitation_raster <- rast("path/to/precipitation.tif")
# 
# # Inspect each raster
# cat("NDVI raster:\n")
# cat("  CRS:", crs(ndvi_raster), "\n")
# cat("  Dimensions:", nrow(ndvi_raster), "x", ncol(ndvi_raster), "\n")
# 
# cat("Elevation raster:\n")
# cat("  CRS:", crs(elevation_raster), "\n")
# cat("  Dimensions:", nrow(elevation_raster), "x", ncol(elevation_raster), "\n")
# 
# # Ensure all rasters have compatible CRS and extent
# if (crs(ndvi_raster) != crs(elevation_raster)) {
#   cat("Warning: Rasters have different CRS\n")
#   # You could reproject here if needed
# }
# 
# # Extract all variables at once
# extracted_multi <- extract_raster_data(
#   raster_files = list(
#     ndvi = ndvi_raster,
#     elevation = elevation_raster,
#     precipitation = precipitation_raster
#   ),
#   study_area = study_area_utm,
#   cell_size = 20000
# )

## ----working-examples, eval = FALSE-------------------------------------------
# # Example 1: WGS84 coordinates (geographic)
# # Use degrees for cell_size
# extracted_wgs84 <- extract_raster_data(
#   raster_files = c(sample = sample_file),
#   study_area = study_area_wgs,  # WGS84 coordinates
#   cell_size = 0.1               # 0.1 degrees (~11km)
# )
# 
# # Example 2: UTM coordinates (projected)
# # Use meters for cell_size
# extracted_utm <- extract_raster_data(
#   raster_files = c(sample = sample_file),
#   study_area = study_area_utm,  # UTM coordinates
#   cell_size = 20000             # 20000 meters (20km)
# )
# 
# # Example 3: No study area (uses raster CRS)
# # Function automatically determines appropriate units
# extracted_auto <- extract_raster_data(
#   raster_files = c(sample = sample_file),
#   cell_size = 0.1               # 0.1 degrees (since raster is WGS84)
# )

## ----configure-topology, eval = TRUE------------------------------------------
# Compute spatial topology for smoothing (default: 2 orders)
topology <- compute_topology(hex_grid)

cat("Topology computed:", length(topology$neighbors[[1]]), "cells\n")
cat("Average distance:", round(topology$avg_distance, 0), "m\n")
cat("Sigma (bandwidth):", round(topology$sigma, 0), "m\n")

cat("Gaussian weights:\n")
cat("Center:", round(topology$weights$center_weight, 3), "\n")
cat("First-order:", round(topology$weights$neighbor_weights[1], 3), "\n") 
cat("Second-order:", round(topology$weights$neighbor_weights[2], 3), "\n")

# Weights applied to individual neighbors, not averages
cat("Weights applied to each neighbor individually\n")

# Neighbor examples
cat("Cell 1:", length(topology$neighbors[[1]][[1]]), "1st,", length(topology$neighbors[[2]][[1]]), "2nd neighbors\n")
cat("Cell 2:", length(topology$neighbors[[1]][[2]]), "1st,", length(topology$neighbors[[2]][[2]]), "2nd neighbors\n")
cat("Cell 3:", length(topology$neighbors[[1]][[3]]), "1st,", length(topology$neighbors[[2]][[3]]), "2nd neighbors\n")

# Note: Visualization code would go here but requires actual grid objects
cat("Neighborhood visualization would show:\n")
cat("- Center cell (red)\n")
cat("- 1st order neighbors (orange)\n") 
cat("- 2nd order neighbors (yellow)\n")

## ----configure-n-order-topology, eval = TRUE----------------------------------
# Compute topology with 3 neighbor orders for more extensive smoothing
topology_3 <- compute_topology(hex_grid, neighbor_orders = 3)

cat("3-Order topology computed:", length(topology_3$neighbors[[1]]), "cells\n")
cat("Neighbor orders:", topology_3$neighbor_orders, "\n")

# Weight computation details
cat("Average distance:", round(topology_3$avg_distance, 0), "m\n")
cat("Sigma (bandwidth):", round(topology_3$sigma, 0), "m\n")

cat("Gaussian weights for 3 orders:\n")
cat("Center:", round(topology_3$weights$center_weight, 3), "\n")
for (i in 1:length(topology_3$weights$neighbor_weights)) {
  cat("Order", i, ":", round(topology_3$weights$neighbor_weights[i], 3), "\n")
}

# Neighbor examples
for (i in 1:min(3, length(topology_3$neighbors[[1]]))) {
  cat("Cell", i, ":")
  for (order in 1:3) {
    n_neighbors <- length(topology_3$neighbors[[order]][[i]])
    cat(" ", n_neighbors, order, "st")
  }
  cat(" neighbors\n")
}

# Custom weights for specific smoothing requirements
topology_custom <- compute_topology(
  hex_grid, 
  neighbor_orders = 4,
  neighbor_weights_param = list(0.5, 0.3, 0.15, 0.05)  # Custom decay pattern
)

cat("Custom 4-order topology with manual weights\n")
cat("Total weight sum:", topology_custom$weights$center_weight + 
    sum(topology_custom$weights$neighbor_weights), "\n")

## ----apply-smoothing, eval = TRUE---------------------------------------------
# Apply spatial smoothing using C++-optimized function (2 orders)
smoothing_results <- smooth_variables(
  variable_values = list(coastal_vegetation = extracted_values),
  neighbors = topology$neighbors,
  weights = topology$weights,
  var_names = c("coastal_vegetation")
)

cat("Spatial smoothing complete\n")
cat("Results:", paste(names(smoothing_results$coastal_vegetation), collapse = ", "), "\n")

# Example smoothing effectiveness analysis:
cat("Smoothing effectiveness:\n")
cat("Original - Mean:", round(mean(extracted_values, na.rm = TRUE), 3), "SD:", round(sd(extracted_values, na.rm = TRUE), 3), "\n")

cat("Smoothed coastal vegetation values:\n")  
cat("  Mean:", round(mean(smoothing_results$coastal_vegetation$weighted_combined, na.rm = TRUE), 3), "\n")
cat("  SD:", round(sd(smoothing_results$coastal_vegetation$weighted_combined, na.rm = TRUE), 3), "\n")

noise_reduction <- (sd(extracted_values, na.rm = TRUE) - 
                   sd(smoothing_results$coastal_vegetation$weighted_combined, na.rm = TRUE)) / 
                   sd(extracted_values, na.rm = TRUE) * 100

cat("Coastal vegetation noise reduction:", round(noise_reduction, 1), "%\n")
cat("(Lower standard deviation = more consistent coastal vegetation patterns)\n")

# The smoothing results contain:
cat("\n=== SMOOTHING RESULTS EXPLAINED ===\n")
cat("- raw: Mean of center cell + all neighbors (unweighted)\n")
cat("- neighbors_1st: Average of first-order neighbors only\n")  
cat("- neighbors_2nd: Average of second-order neighbors only\n")
cat("- weighted_combined: Weighted average (center + all neighbors) ★\n")
cat("\n★ weighted_combined is typically the best result for noise reduction\n")

## ----apply-n-order-smoothing, eval = TRUE-------------------------------------
# Apply 3-order smoothing using the new N-order system
smoothing_results_3 <- smooth_variables(
  variable_values = list(coastal_vegetation = extracted_values),
  neighbors = topology_3$neighbors,
  weights = topology_3$weights,
  var_names = c("coastal_vegetation")
)

cat("3-Order spatial smoothing complete\n")
cat("Results:", paste(names(smoothing_results_3$coastal_vegetation), collapse = ", "), "\n")

# Compare 2-order vs 3-order smoothing
cat("\n=== COMPARING SMOOTHING APPROACHES ===\n")

# 2-order results
cat("2-Order smoothing:\n")
cat("  Mean:", round(mean(smoothing_results$coastal_vegetation$weighted_combined, na.rm = TRUE), 3), "\n")
cat("  SD:", round(sd(smoothing_results$coastal_vegetation$weighted_combined, na.rm = TRUE), 3), "\n")

# 3-order results
cat("3-Order smoothing:\n")
cat("  Mean:", round(mean(smoothing_results_3$coastal_vegetation$weighted_combined, na.rm = TRUE), 3), "\n")
cat("  SD:", round(sd(smoothing_results_3$coastal_vegetation$weighted_combined, na.rm = TRUE), 3), "\n")

# Additional 3-order results available
cat("\n3-Order specific results:\n")
cat("- neighbors_1st:", round(mean(smoothing_results_3$coastal_vegetation$neighbors_1st, na.rm = TRUE), 3), "\n")
cat("- neighbors_2nd:", round(mean(smoothing_results_3$coastal_vegetation$neighbors_2nd, na.rm = TRUE), 3), "\n")
cat("- neighbors_3rd:", round(mean(smoothing_results_3$coastal_vegetation$neighbors_3rd, na.rm = TRUE), 3), "\n")

# Calculate additional noise reduction from 3-order smoothing
noise_reduction_3 <- (sd(extracted_values, na.rm = TRUE) - 
                     sd(smoothing_results_3$coastal_vegetation$weighted_combined, na.rm = TRUE)) / 
                     sd(extracted_values, na.rm = TRUE) * 100

cat("3-Order noise reduction:", round(noise_reduction_3, 1), "%\n")
cat("Additional reduction over 2-order:", round(noise_reduction_3 - noise_reduction, 1), "%\n")

## ----combine-results, eval = TRUE---------------------------------------------
# Combine the hexagonal grid with smoothing results
hex_grid_with_results <- hex_grid
hex_grid_with_results$coastal_vegetation_raw <- smoothing_results$coastal_vegetation$raw
hex_grid_with_results$coastal_vegetation_neighbors_1st <- smoothing_results$coastal_vegetation$neighbors_1st
hex_grid_with_results$coastal_vegetation_neighbors_2nd <- smoothing_results$coastal_vegetation$neighbors_2nd
hex_grid_with_results$coastal_vegetation_weighted_combined <- smoothing_results$coastal_vegetation$weighted_combined

# Transform to WGS84 for plotting
hex_grid_with_results_wgs84 <- st_transform(hex_grid_with_results, crs = 4326)

cat("Results combined with grid successfully\n")
cat("Available columns:", paste(names(hex_grid_with_results), collapse = ", "), "\n")

# Analyze the smoothing effect
cat("\n=== SMOOTHING ANALYSIS ===\n")
cat("Understanding the smoothing results:\n")
cat("- raw: Mean of center cell + all neighbors (unweighted)\n")
cat("- neighbors_1st: Mean of first-order neighbors only\n")
cat("- neighbors_2nd: Mean of second-order neighbors only\n")
cat("- weighted_combined: Weighted average (center + all neighbors, weights applied to each neighbor)\n\n")

# Example analysis output:
cat("Original values - Count:", length(extracted_values), "Mean:", round(mean(extracted_values, na.rm = TRUE), 4), "SD:", round(sd(extracted_values, na.rm = TRUE), 4), "\n")
cat("Smoothed values - Count:", length(smoothing_results$coastal_vegetation$weighted_combined), "Mean:", round(mean(smoothing_results$coastal_vegetation$weighted_combined, na.rm = TRUE), 4), "SD:", round(sd(smoothing_results$coastal_vegetation$weighted_combined, na.rm = TRUE), 4), "\n")
cat("Variance reduction:", round(noise_reduction, 1), "%\n")

## ----visualize-results, eval = TRUE-------------------------------------------
# Create spectacular before/after visualization showing the power of spatial smoothing
cat("=== CREATING BEFORE/AFTER VISUALIZATION ===\n")

# Create visualization plots
par(mfrow = c(2, 2))

# Plot 1: Original extracted values
plot(st_geometry(hex_grid_wgs84), 
     col = plot_colors[cut(extracted_values, 100)],
     main = "Original Extracted Values",
     border = "white", lwd = 0.01)

# Plot 2: Smoothed values (1st order)
plot(st_geometry(hex_grid_wgs84), 
     col = plot_colors[cut(smoothing_results$coastal_vegetation$neighbors_1st, 100)],
     main = "Smoothed (1st Order Neighbors)",
     border = "white", lwd = 0.01)

# Plot 3: Smoothed values (2nd order)
plot(st_geometry(hex_grid_wgs84), 
     col = plot_colors[cut(smoothing_results$coastal_vegetation$neighbors_2nd, 100)],
     main = "Smoothed (2nd Order Neighbors)",
     border = "white", lwd = 0.01)

# Plot 4: Final weighted combined result
plot(st_geometry(hex_grid_wgs84), 
     col = plot_colors[cut(smoothing_results$coastal_vegetation$weighted_combined, 100)],
     main = "Final Weighted Combined Result",
     border = "white", lwd = 0.01)

cat("All visualization plots created successfully!\n")
cat("Progression shows:\n")
cat("- Top left: Raw extracted data\n")
cat("- Top right: 1st order neighbor smoothing\n")
cat("- Bottom left: 2nd order neighbor smoothing\n")
cat("- Bottom right: Final weighted combined result\n")

## ----multiple-variables, eval = FALSE-----------------------------------------
# # You can work with multiple raster files and variables
# raster_files <- c(
#   ndvi = "path/to/your/ndvi.tif",
#   elevation = "path/to/your/elevation.tif",
#   precipitation = "path/to/your/precipitation.tif"
# )
# 
# # Extract multiple variables
# extracted_multi <- extract_raster_data(
#   raster_files = raster_files,
#   study_area = study_area_utm,
#   cell_size = cell_size
# )
# 
# # Apply smoothing to all variables
# variable_values <- list(
#   ndvi = extracted_multi$data$ndvi,
#   elevation = extracted_multi$data$elevation,
#   precipitation = extracted_multi$data$precipitation
# )
# 
# smoothing_multi <- smooth_variables(
#   variable_values = variable_values,
#   neighbors = topology$neighbors,
#   weights = topology$weights,
#   var_names = c("ndvi", "elevation", "precipitation")
# )
# 
# # Each variable gets the same smoothing treatment
# cat("Variables processed:", paste(names(smoothing_multi), collapse = ", "), "\n")

## ----hexagon-utilities, eval = TRUE-------------------------------------------
# The package provides helper functions for hexagon measurements
# IMPORTANT: cell_size_flat is in METERS (UTM coordinates)
cell_size_flat <- 20000  # 20km flat-to-flat distance in meters

# Convert between different hexagon measurements
edge_length <- hex_flat_to_edge(cell_size_flat)
circumradius <- hex_flat_to_circumradius(cell_size_flat)

cat("Hexagon measurements for", cell_size_flat, "m flat-to-flat distance:\n")
cat("     - Edge length:", round(edge_length, 2), "m\n")
cat("     - Circumradius:", round(circumradius, 2), "m\n")

# Convert back
flat_distance_from_edge <- hex_edge_to_flat(edge_length)
flat_distance_from_circumradius <- hex_circumradius_to_flat(circumradius)

cat("Verification:\n")
cat("Flat distance from edge:", round(flat_distance_from_edge, 2), "m\n")
cat("Flat distance from circumradius:", round(flat_distance_from_circumradius, 2), "m\n")

