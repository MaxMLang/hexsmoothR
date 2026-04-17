# hexsmoothR 0.2.0

## New features

* `create_grid()` and `compute_topology()` now auto-detect a sensible UTM
  projection via `get_utm_crs()` when `projection_crs = NULL`.
* `find_hex_cell_size_for_target_cells()` derives sensible
  `cell_size_min`/`cell_size_max` defaults from the study-area bounding box
  when those arguments are `NULL`.
* New global option `hexsmoothR.verbose` (default `TRUE`) controls all
  informational console output via the internal `hex_msg()` helper.

## Bug fixes

* `extract_raster_data()` now reprojects the hex grid to each raster's CRS
  before calling `exactextractr::exact_extract()`, fixing silent CRS
  mismatches with multi-raster inputs. The cropped study area is also
  rebuilt per raster.
* `extract_raster_data()` returns the (possibly sampled) hex grid in the
  result list as `hex_grid`, and the "Sampled X from Y" log line is now
  correct.
* `get_utm_crs()` correctly handles longitudes near 180 and auto-projects
  inputs to WGS84 before computing the centroid.
* `smooth_variables()`: the `raw` element of the output now stores the
  original (unsmoothed) value of the centre cell, as documented.
* C++ `Rcpp::stop()` messages no longer print literal `%d` placeholders.
* `compute_topology()` order-N BFS rewritten using set operations; large
  grids and higher orders are noticeably faster.

## Documentation

* Vignette `hexsmoothR-complete-guide` has been rewritten end-to-end and is
  now driven by the bundled `inst/extdata/default.tif` sample raster.
* Added `NEWS.md`, `BugReports`, and a second `URL` for the GitHub repo.

## Internal

* Test suite reorganised into focused files
  (`test-create-grid.R`, `test-compute-topology.R`,
  `test-smooth-variables.R`, `test-utils.R`, `test-get-utm-crs.R`,
  `test-vignette-runs.R`, `test-vignette-workflow.R`,
  `test-n-order-smoothing.R`).
* `useDynLib(hexsmoothR, .registration = TRUE)` added so registered C++
  symbols are correctly resolved.

# hexsmoothR 0.1.0

* Initial release.
