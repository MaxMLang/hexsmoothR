# Package index

## Grids

Build and tune hexagonal / rectangular grids over a study area.

- [`create_grid()`](https://maxmlang.github.io/hexsmoothR/reference/create_grid.md)
  : Create hexagonal or square grids for spatial analysis
- [`create_single_grid()`](https://maxmlang.github.io/hexsmoothR/reference/create_single_grid.md)
  : Build the grid for a single study area (internal)
- [`find_hex_cell_size_for_target_cells()`](https://maxmlang.github.io/hexsmoothR/reference/find_hex_cell_size_for_target_cells.md)
  : Find a hexagonal cell size that yields approximately a target number
  of cells
- [`get_utm_crs()`](https://maxmlang.github.io/hexsmoothR/reference/get_utm_crs.md)
  : Get an appropriate UTM CRS for a study area

## Raster extraction

Pull raster values into grid cells.

- [`extract_raster_data()`](https://maxmlang.github.io/hexsmoothR/reference/extract_raster_data.md)
  : Extract raster data into a hexagonal grid

## Topology and smoothing

Compute neighbour topology and apply N-order spatial smoothing.

- [`compute_topology()`](https://maxmlang.github.io/hexsmoothR/reference/compute_topology.md)
  : Compute spatial topology and weights for hexagonal smoothing
- [`smooth_variables()`](https://maxmlang.github.io/hexsmoothR/reference/smooth_variables.md)
  : Apply spatial smoothing to variables on a hexagonal grid

## Hex measurement helpers

- [`hex_flat_to_edge()`](https://maxmlang.github.io/hexsmoothR/reference/hex_flat_to_edge.md)
  [`hex_flat_to_circumradius()`](https://maxmlang.github.io/hexsmoothR/reference/hex_flat_to_edge.md)
  [`hex_edge_to_flat()`](https://maxmlang.github.io/hexsmoothR/reference/hex_flat_to_edge.md)
  [`hex_circumradius_to_flat()`](https://maxmlang.github.io/hexsmoothR/reference/hex_flat_to_edge.md)
  : Convert between hexagon measurements

## Internal

- [`hex_msg()`](https://maxmlang.github.io/hexsmoothR/reference/hex_msg.md)
  : Internal verbose-print helper
- [`ordinal_suffix()`](https://maxmlang.github.io/hexsmoothR/reference/ordinal_suffix.md)
  : Ordinal suffix for an integer (internal)
- [`validate_smoothing_inputs()`](https://maxmlang.github.io/hexsmoothR/reference/validate_smoothing_inputs.md)
  : Validate inputs to the smoothing functions (internal)
- [`smooth_variables_r_fallback()`](https://maxmlang.github.io/hexsmoothR/reference/smooth_variables_r_fallback.md)
  : R fallback for N-order spatial smoothing
- [`estimate_avg_neighbor_distance()`](https://maxmlang.github.io/hexsmoothR/reference/estimate_avg_neighbor_distance.md)
  : Estimate average distance between first-order neighbours (internal)
- [`compute_neighbors_n_orders()`](https://maxmlang.github.io/hexsmoothR/reference/compute_neighbors_n_orders.md)
  : Compute neighbours for N orders (internal, vectorised set-based BFS)
- [`compute_weights_n_orders()`](https://maxmlang.github.io/hexsmoothR/reference/compute_weights_n_orders.md)
  : Compute weights for N orders (internal)
