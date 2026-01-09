# Find optimal hexagonal cell size for target number of cells

Uses binary search to find the hexagonal cell size (flat-to-flat
distance) that produces approximately the target number of hexagons.

## Usage

``` r
find_hex_cell_size_for_target_cells(
  study_area,
  target_cells,
  cell_size_min = 0.001,
  cell_size_max = 10,
  tol = 0.05,
  max_iter = 20
)
```

## Arguments

- study_area:

  sf object (polygons)

- target_cells:

  target number of hexagons

- cell_size_min:

  minimum cell size to try (in map units)

- cell_size_max:

  maximum cell size to try (in map units)

- tol:

  tolerance for convergence (fraction of target_cells)

- max_iter:

  maximum iterations for binary search

## Value

cell size (flat-to-flat distance) that gives closest to target_cells
hexagons
