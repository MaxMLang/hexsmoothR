# Smoke test: the vignette code itself must run end-to-end without error.
# We re-implement the vignette as a single script so a CI run catches drift
# even if knitr / rmarkdown aren't installed.

library(testthat)
library(hexsmoothR)
library(sf)
library(terra)

test_that("vignette workflow runs end-to-end", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  options(hexsmoothR.verbose = FALSE)
  set.seed(1)

  ndvi_path <- system.file("extdata", "default.tif", package = "hexsmoothR")
  expect_true(file.exists(ndvi_path))

  ndvi <- rast(ndvi_path)
  expect_s4_class(ndvi, "SpatRaster")

  study_area_wgs <- st_as_sf(st_as_sfc(st_bbox(ndvi)))
  st_crs(study_area_wgs) <- 4326

  utm_crs <- get_utm_crs(study_area_wgs)
  expect_match(utm_crs, "^EPSG:32[67][0-9]{2}$")

  study_area_utm <- st_transform(study_area_wgs, utm_crs)

  hex_grid <- create_grid(
    study_area = study_area_utm,
    cell_size  = 25000,
    type       = "hexagonal"
  )
  expect_s3_class(hex_grid, "sf")
  expect_gt(nrow(hex_grid), 10)

  target_size <- find_hex_cell_size_for_target_cells(
    study_area    = study_area_utm,
    target_cells  = 250,
    cell_size_min = 5000,
    cell_size_max = 100000
  )
  expect_gt(target_size, 5000)
  expect_lt(target_size, 100000)

  extracted <- extract_raster_data(
    raster_files = list(ndvi = ndvi),
    hex_grid     = hex_grid
  )
  expect_named(extracted, c("data", "hex_grid", "cell_size", "extent",
                            "variables", "n_cells"))
  expect_equal(extracted$variables, "ndvi")
  expect_s3_class(extracted$hex_grid, "sf")
  expect_equal(nrow(extracted$data), nrow(hex_grid))
  expect_true(any(!is.na(extracted$data$ndvi)))

  topology <- compute_topology(hex_grid, neighbor_orders = 2)
  expect_equal(topology$neighbor_orders, 2)
  expect_equal(
    topology$weights$center_weight + sum(topology$weights$neighbor_weights),
    1, tolerance = 1e-10
  )

  topo_custom <- compute_topology(
    hex_grid,
    neighbor_orders        = 3,
    neighbor_weights_param = list(0.5, 0.3, 0.1)
  )
  expect_length(topo_custom$weights$neighbor_weights, 3)

  result <- smooth_variables(
    variable_values = list(ndvi = extracted$data$ndvi),
    neighbors       = topology$neighbors,
    weights         = topology$weights
  )
  expect_named(result, "ndvi")
  expect_named(
    result$ndvi,
    c("raw", "weighted_combined", "neighbors_1st", "neighbors_2nd"),
    ignore.order = TRUE
  )
  expect_equal(result$ndvi$raw, extracted$data$ndvi)
  expect_false(isTRUE(all.equal(
    result$ndvi$raw, result$ndvi$weighted_combined
  )))

  sds <- vapply(2:4, function(k) {
    topo_k <- compute_topology(hex_grid, neighbor_orders = k)
    out <- smooth_variables(
      variable_values = list(ndvi = extracted$data$ndvi),
      neighbors       = topo_k$neighbors,
      weights         = topo_k$weights
    )
    sd(out$ndvi$weighted_combined, na.rm = TRUE)
  }, numeric(1))
  expect_true(all(is.finite(sds)))

  multi <- smooth_variables(
    variable_values = list(
      ndvi         = extracted$data$ndvi,
      ndvi_squared = extracted$data$ndvi^2
    ),
    neighbors = topology$neighbors,
    weights   = topology$weights
  )
  expect_named(multi, c("ndvi", "ndvi_squared"))
})
