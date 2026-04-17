# Extract raster data into a hexagonal grid

Primary function for extracting raster values into hexagonal grid cells.
CRS transformations are handled automatically: each raster is
reprojected (or the grid is reprojected to the raster's CRS) so that the
underlying call to \[\`exactextractr::exact_extract()\`\] always sees
matching CRSs.

## Usage

``` r
extract_raster_data(
  raster_files,
  study_area = NULL,
  cell_size = NULL,
  hex_grid = NULL,
  sample_fraction = 1,
  random_seed = 42,
  fun = "mean"
)
```

## Arguments

- raster_files:

  Named character vector of file paths OR named list of
  \`terra::SpatRaster\` objects.

- study_area:

  Optional sf polygon used for cropping each raster and to define the
  grid CRS.

- cell_size:

  Hex cell size, in the units of the grid CRS. Required when
  \`hex_grid\` is not supplied.

- hex_grid:

  Optional sf hexagonal grid to use instead of creating one.

- sample_fraction:

  Fraction of grid cells to keep (default 1).

- random_seed:

  Seed for reproducible sampling.

- fun:

  Aggregation function passed to \`exactextractr::exact_extract()\`
  (default \`"mean"\`).

## Value

List with components

- \`data\`:

  Data frame with \`cell_id\`, \`x\`, \`y\` and one column per raster.

- \`hex_grid\`:

  The sf grid that was used (sampled, if applicable).

- \`cell_size\`:

  The cell size used.

- \`extent\`:

  Extent of the first raster (after cropping).

- \`variables\`:

  Names of the rasters.

- \`n_cells\`:

  Number of cells in \`data\`.

## Input flexibility

\`raster_files\` may be either a named character vector of file paths or
a named list of \`terra::SpatRaster\` objects. All inputs may be in
different CRSs from one another and from the grid - the function handles
cropping and transformation per raster.

## CRS handling

\- If \`study_area\` is supplied, the grid is created in the study
area's CRS. - Otherwise the grid is created in the first raster's CRS. -
For each raster, the grid is transformed to the raster's CRS before
extraction (so cell-size units are honoured exactly once, in the grid
CRS).
