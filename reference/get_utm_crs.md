# Get appropriate UTM CRS for a study area

Automatically determines the appropriate UTM coordinate reference system
(CRS) for a given study area based on its centroid location.

## Usage

``` r
get_utm_crs(study_area)
```

## Arguments

- study_area:

  sf object representing the study area

## Value

Character string with the UTM CRS (e.g., "EPSG:32630" for UTM zone 30N)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)

# Create a simple study area
study_area <- st_sf(geometry = st_sfc(
  st_polygon(list(matrix(c(-5, 35, 5, 35, 5, 45, -5, 45, -5, 35), ncol = 2, byrow = TRUE))),
  crs = 4326
))

# Get appropriate UTM CRS
utm_crs <- get_utm_crs(study_area)
print(utm_crs)  # Should return something like "EPSG:32630"
} # }
```
