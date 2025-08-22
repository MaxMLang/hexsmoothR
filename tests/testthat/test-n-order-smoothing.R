# Test N-order smoothing functionality
# This test file verifies that the new configurable neighbor orders work correctly
# while maintaining full backward compatibility with existing code.

library(testthat)
library(hexsmoothR)
library(sf)

# Create a simple test grid
create_test_grid <- function() {
  # Create a small study area in WGS84
  study_area_wgs <- st_sf(geometry = st_sfc(
    st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE))),
    crs = 4326
  ))
  
  # Transform to UTM (like the vignette workflow)
  utm_crs <- get_utm_crs(study_area_wgs)
  study_area_utm <- st_transform(study_area_wgs, utm_crs)
  
  # Create a hexagonal grid with the same cell size as vignette (20km)
  grid <- create_grid(study_area_utm, cell_size = 20000, type = "hexagonal")
  return(grid)
}

test_that("compute_topology works with 2 neighbor orders", {
  grid <- create_test_grid()
  
  # Test default behavior (2 orders)
  topology <- compute_topology(grid)
  
  # Verify new fields exist
  expect_true("neighbors" %in% names(topology))
  expect_true("weights" %in% names(topology))
  expect_true("neighbor_orders" %in% names(topology))
  
  # Verify neighbor_orders is 2
  expect_equal(topology$neighbor_orders, 2)
  
  # Verify neighbors list has 2 elements
  expect_equal(length(topology$neighbors), 2)
  
  # Verify weights structure
  expect_true("center_weight" %in% names(topology$weights))
  expect_true("neighbor_weights" %in% names(topology$weights))
  expect_equal(length(topology$weights$neighbor_weights), 2)
})

test_that("compute_topology works with 3 neighbor orders", {
  grid <- create_test_grid()
  
  # Test 3 orders
  topology <- compute_topology(grid, neighbor_orders = 3)
  
  # Verify structure
  expect_equal(topology$neighbor_orders, 3)
  expect_equal(length(topology$neighbors), 3)
  expect_equal(length(topology$weights$neighbor_weights), 3)
  
  # Verify no backward compatibility fields for 3 orders
  expect_false("first_order_neighbors" %in% names(topology))
  expect_false("second_order_neighbors" %in% names(topology))
  
  # Verify neighbors structure
  expect_true(is.list(topology$neighbors))
  expect_equal(length(topology$neighbors), 3)
  
  # Verify weights decay appropriately
  weights <- topology$weights$neighbor_weights
  expect_true(weights[1] > weights[2])  # First order > second order
  expect_true(weights[2] > weights[3])  # Second order > third order
})

test_that("compute_topology works with custom weights", {
  grid <- create_test_grid()
  
  # Test custom weights for 3 orders
  custom_weights <- list(0.6, 0.3, 0.1)
  topology <- compute_topology(grid, neighbor_orders = 3, neighbor_weights_param = custom_weights)
  
  # Verify custom weights are used (normalized)
  # Original weights: [0.6, 0.3, 0.1] with center_weight = 1.0, total = 2.0
  # After normalization: each weight / 2.0
  expect_equal(topology$weights$neighbor_weights, c(0.30, 0.15, 0.05), tolerance = 0.001)
  
  # Verify weights are normalized
  total_weight <- topology$weights$center_weight + sum(topology$weights$neighbor_weights)
  expect_equal(total_weight, 1.0, tolerance = 1e-10)
})

test_that("smooth_variables works with 2 neighbor orders", {
  grid <- create_test_grid()
  topology <- compute_topology(grid)  # 2 orders
  
  # Create test data
  variable_values <- list(
    test_var = runif(nrow(grid), min = 1.0, max = 10.0)
  )
  
  # Test N-order call with 2 orders
  results <- smooth_variables(
    variable_values = variable_values,
    neighbors = topology$neighbors,
    weights = topology$weights,
    var_names = c("test_var")
  )
  
  # Verify results structure
  expect_true(is.numeric(results$test_var$raw))
  expect_true(is.numeric(results$test_var$neighbors_1st))
  expect_true(is.numeric(results$test_var$neighbors_2nd))
  expect_true(is.numeric(results$test_var$weighted_combined))
  
  # Verify results have the expected length
  expect_equal(length(results$test_var$raw), nrow(grid))
  expect_equal(length(results$test_var$neighbors_1st), nrow(grid))
  expect_equal(length(results$test_var$neighbors_2nd), nrow(grid))
})

test_that("smooth_variables works with 3 neighbor orders", {
  grid <- create_test_grid()
  topology <- compute_topology(grid, neighbor_orders = 3)
  
  # Create test data
  variable_values <- list(
    test_var = runif(nrow(grid), min = 1.0, max = 10.0)
  )
  
  # Test 3-order smoothing
  results <- smooth_variables(
    variable_values = variable_values,
    neighbors = topology$neighbors,
    weights = topology$weights,
    var_names = c("test_var")
  )
  
  # Verify structure
  expect_true("neighbors_1st" %in% names(results$test_var))
  expect_true("neighbors_2nd" %in% names(results$test_var))
  expect_true("neighbors_3rd" %in% names(results$test_var))
  expect_true("weighted_combined" %in% names(results$test_var))
  
  # Verify results are numeric
  expect_true(is.numeric(results$test_var$neighbors_1st))
  expect_true(is.numeric(results$test_var$neighbors_2nd))
  expect_true(is.numeric(results$test_var$neighbors_3rd))
  expect_true(is.numeric(results$test_var$weighted_combined))
})

test_that("parameter validation works correctly", {
  grid <- create_test_grid()
  
  # Test invalid neighbor_orders
  expect_error(compute_topology(grid, neighbor_orders = 0))
  expect_error(compute_topology(grid, neighbor_orders = -1))
  expect_error(compute_topology(grid, neighbor_orders = 1.5))
  
  # Test invalid neighbor_weights_param
  expect_error(compute_topology(grid, neighbor_orders = 3, 
                               neighbor_weights_param = list(0.5, 0.3)))  # Wrong length
  expect_error(compute_topology(grid, neighbor_orders = 2, 
                               neighbor_weights_param = list("invalid", 0.3)))  # Non-numeric
})

test_that("memory and performance scaling works reasonably", {
  grid <- create_test_grid()
  
  # Test different orders and measure performance
  orders_to_test <- c(2, 3, 4)
  
  for (order in orders_to_test) {
    start_time <- Sys.time()
    topology <- compute_topology(grid, neighbor_orders = order)
    end_time <- Sys.time()
    
    # Verify structure
    expect_equal(topology$neighbor_orders, order)
    expect_equal(length(topology$neighbors), order)
    expect_equal(length(topology$weights$neighbor_weights), order)
    
    # Verify computation time is reasonable (should not be exponential)
    computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    expect_true(computation_time < 10)  # Should complete within 10 seconds
    
    cat(sprintf("Order %d completed in %.3f seconds\n", order, computation_time))
  }
})

test_that("weight normalization works correctly", {
  grid <- create_test_grid()
  
  # Test with custom weights that don't sum to 1
  custom_weights <- list(0.8, 0.4, 0.2)  # Sum = 1.4
  topology <- compute_topology(grid, neighbor_orders = 3, neighbor_weights_param = custom_weights)
  
  # Verify weights are normalized
  total_weight <- topology$weights$center_weight + sum(topology$weights$neighbor_weights)
  expect_equal(total_weight, 1.0, tolerance = 1e-10)
  
  # Verify relative proportions are maintained
  original_ratio <- custom_weights[[1]] / custom_weights[[2]]
  normalized_ratio <- topology$weights$neighbor_weights[1] / topology$weights$neighbor_weights[2]
  expect_equal(original_ratio, normalized_ratio, tolerance = 1e-10)
})

test_that("edge cases are handled correctly", {
  grid <- create_test_grid()
  
  # Test with very high orders (should not crash)
  expect_no_error(compute_topology(grid, neighbor_orders = 10))
  
  # Test with single cell grid (edge case)
  single_cell_grid <- grid[1, ]
  expect_no_error(compute_topology(single_cell_grid, neighbor_orders = 3))
  
  # Test with empty grid (should handle gracefully)
  empty_grid <- grid[integer(0), ]
  result <- suppressWarnings(compute_topology(empty_grid, neighbor_orders = 2))
  expect_equal(length(result), 0)  # Should return empty list
}) 