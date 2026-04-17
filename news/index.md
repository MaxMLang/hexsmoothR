# Changelog

## hexsmoothR 0.2.0

### New features

- [`create_grid()`](https://maxmlang.github.io/hexsmoothR/reference/create_grid.md)
  and
  [`compute_topology()`](https://maxmlang.github.io/hexsmoothR/reference/compute_topology.md)
  now auto-detect a sensible UTM projection via
  [`get_utm_crs()`](https://maxmlang.github.io/hexsmoothR/reference/get_utm_crs.md)
  when `projection_crs = NULL`.
- [`find_hex_cell_size_for_target_cells()`](https://maxmlang.github.io/hexsmoothR/reference/find_hex_cell_size_for_target_cells.md)
  derives sensible `cell_size_min`/`cell_size_max` defaults from the
  study-area bounding box when those arguments are `NULL`.
- New global option `hexsmoothR.verbose` (default `TRUE`) controls all
  informational console output via the internal
  [`hex_msg()`](https://maxmlang.github.io/hexsmoothR/reference/hex_msg.md)
  helper.

### Bug fixes

- [`extract_raster_data()`](https://maxmlang.github.io/hexsmoothR/reference/extract_raster_data.md)
  now reprojects the hex grid to each raster’s CRS before calling
  [`exactextractr::exact_extract()`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html),
  fixing silent CRS mismatches with multi-raster inputs. The cropped
  study area is also rebuilt per raster.
- [`extract_raster_data()`](https://maxmlang.github.io/hexsmoothR/reference/extract_raster_data.md)
  returns the (possibly sampled) hex grid in the result list as
  `hex_grid`, and the “Sampled X from Y” log line is now correct.
- [`get_utm_crs()`](https://maxmlang.github.io/hexsmoothR/reference/get_utm_crs.md)
  correctly handles longitudes near 180 and auto-projects inputs to
  WGS84 before computing the centroid.
- [`smooth_variables()`](https://maxmlang.github.io/hexsmoothR/reference/smooth_variables.md):
  the `raw` element of the output now stores the original (unsmoothed)
  value of the centre cell, as documented.
- C++ `Rcpp::stop()` messages no longer print literal `%d` placeholders.
- [`compute_topology()`](https://maxmlang.github.io/hexsmoothR/reference/compute_topology.md)
  order-N BFS rewritten using set operations; large grids and higher
  orders are noticeably faster.

### Documentation

- Vignette `hexsmoothR-complete-guide` has been rewritten end-to-end and
  is now driven by the bundled `inst/extdata/default.tif` sample raster.
- Added `NEWS.md`, `BugReports`, and a second `URL` for the GitHub repo.

### Internal

- Test suite reorganised into focused files (`test-create-grid.R`,
  `test-compute-topology.R`, `test-smooth-variables.R`, `test-utils.R`,
  `test-get-utm-crs.R`, `test-vignette-runs.R`,
  `test-vignette-workflow.R`, `test-n-order-smoothing.R`).
- `useDynLib(hexsmoothR, .registration = TRUE)` added so registered C++
  symbols are correctly resolved.

## hexsmoothR 0.1.0

- Initial release.
