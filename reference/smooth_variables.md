# Apply spatial smoothing to variables using hexagonal grid topology

This function applies spatial smoothing to variables using a hexagonal
grid topology. It can be used independently after data extraction, or as
part of the complete workflow. The function automatically uses the C++
implementation when available and falls back to R if needed.

## Usage

``` r
smooth_variables(
  variable_values,
  neighbors,
  weights,
  hex_indices = NULL,
  var_names
)
```

## Arguments

- variable_values:

  List of numeric vectors containing variable values for each grid cell

- neighbors:

  List of neighbor indices for each order (from compute_topology)

- weights:

  List containing center_weight and neighbor_weights for each order

- hex_indices:

  Vector of hexagon indices to process (default: all cells)

- var_names:

  Character vector of variable names

## Value

List containing smoothing results for each variable with the following
components:

- raw: Weighted average of center cell and all neighbors

- neighbors_Nst: Mean of neighbors at order N

- weighted_combined: Weighted average of center cell and all neighbors

## Examples

``` r
# After creating a grid and computing topology
# grid_sf <- create_grid(study_area, cell_size = 10000)
# topology <- compute_topology(grid_sf, neighbor_orders = 3)  # 3 orders

# Extract your variables (example)
# variable_values <- list(
#   ndvi = c(0.5, 0.6, 0.4, 0.7, 0.3),
#   elevation = c(100, 120, 90, 140, 80)
# )

# Apply smoothing with N orders
# smoothing_results <- smooth_variables(
#   variable_values = variable_values,
#   neighbors = topology$neighbors,
#   weights = topology$weights,
#   var_names = c("ndvi", "elevation")
# )
```
