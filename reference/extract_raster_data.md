# Extract raster data for a hex grid using exactextractr

This is the PRIMARY function for extracting raster data to hexagonal
grids. It automatically handles CRS transformations and creates
consistent hexagonal grids.

## Usage

``` r
extract_raster_data(
  raster_files,
  study_area = NULL,
  cell_size = NULL,
  hex_grid = NULL,
  sample_fraction = 1,
  random_seed = 42
)
```

## Arguments

- raster_files:

  named vector of raster file paths OR named list of terra::rast objects

- study_area:

  sf polygon (optional, for cropping and grid CRS reference)

- cell_size:

  hex cell size (in units of the target CRS - see details)

- hex_grid:

  optional sf hexagonal grid to use instead of creating a new one

- sample_fraction:

  fraction of raster cells to use (default 1)

- random_seed:

  for reproducibility

## Value

list: data frame with cell_id, x, y, and values for each raster

## Details

\*\*CRS Handling:\*\* - If \`study_area\` is provided, the function
creates a grid in the study area's CRS - If no \`study_area\` is
provided, the function creates a grid in the raster's CRS - The function
automatically crops rasters to the study area (if provided) - All
coordinate transformations are handled transparently

\*\*Grid Creation:\*\* - Creates hexagonal grids using
\`sf::st_make_grid\` for robustness - Automatically adjusts cell sizes
for geographic coordinates (degrees) - Can use a pre-existing grid via
the \`hex_grid\` parameter

\*\*Input Flexibility:\*\* - Accepts both file paths (character strings)
and \`terra::rast\` objects - When passing \`terra::rast\` objects, you
can inspect CRS, dimensions, and data before extraction - Useful for
working with in-memory rasters or pre-processed data

\*\*CRS Handling Details:\*\*

The function intelligently handles coordinate reference systems:

1\. \*\*Study Area Provided\*\*: Grid is created in the study area's
CRS - Example: If study_area is in UTM (EPSG:32630), grid uses UTM
coordinates - Cell sizes should be in meters (e.g., 20000 for 20km
hexagons)

2\. \*\*No Study Area\*\*: Grid is created in the raster's CRS -
Example: If raster is in WGS84 (EPSG:4326), grid uses geographic
coordinates - Cell sizes should be in degrees (e.g., 0.1 for ~11km
hexagons)

3\. \*\*Automatic Adjustments\*\*: - Large cell sizes (\>1) in
geographic coordinates trigger warnings - Cell sizes are automatically
adjusted for geographic coordinates if needed

\*\*Input Types:\*\*

The function accepts two types of input for \`raster_files\`:

1\. \*\*File Paths\*\*: Named character vector of file paths “\`r
raster_files = c(ndvi = "path/to/ndvi.tif", elevation =
"path/to/elevation.tif") “\`

2\. \*\*Terra Raster Objects\*\*: Named list of terra::rast objects “\`r
\# Load rasters first to inspect them ndvi_rast \<-
rast("path/to/ndvi.tif") elev_rast \<- rast("path/to/elevation.tif")

\# Check CRS, dimensions, etc. print(crs(ndvi_rast))
print(dim(ndvi_rast))

\# Pass to function raster_files = list(ndvi = ndvi_rast, elevation =
elev_rast) “\`

\*\*Workflow Recommendation:\*\*

For consistent results, always provide a \`study_area\` in your desired
CRS: “\`r \# Create study area in UTM for accurate measurements
study_area_utm \<- st_transform(study_area_wgs, utm_crs)

\# Extract with 20km hexagons (in UTM coordinates) extracted \<-
extract_raster_data( raster_files = c(var1 = "path/to/raster.tif"),
study_area = study_area_utm, cell_size = 20000 \# 20km in meters )

\# Then apply smoothing smoothed \<- smooth_variables(...) “\`
