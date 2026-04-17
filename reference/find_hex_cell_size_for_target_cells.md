# Find a hexagonal cell size that yields approximately a target number of cells

Uses a binary search over \`cell_size_min\` / \`cell_size_max\` to find
the flat-to-flat distance that produces approximately \`target_cells\`
hexagons.

## Usage

``` r
find_hex_cell_size_for_target_cells(
  study_area,
  target_cells,
  cell_size_min = NULL,
  cell_size_max = NULL,
  tol = 0.05,
  max_iter = 20,
  projection_crs = NULL
)
```

## Arguments

- study_area:

  sf object (polygons).

- target_cells:

  Desired number of hexagons (positive integer).

- cell_size_min:

  Minimum cell size to try. Default \`NULL\`, in which case the search
  range is derived from the bounding box of \`study_area\` so it works
  for either projected (metres) or geographic (degrees) input.

- cell_size_max:

  Maximum cell size to try. See \`cell_size_min\`.

- tol:

  Convergence tolerance (fraction of \`target_cells\`).

- max_iter:

  Maximum binary-search iterations.

- projection_crs:

  Optional CRS for grid construction (passed to \[\`create_grid()\`\]).

## Value

Cell size (flat-to-flat distance) closest to \`target_cells\`.
