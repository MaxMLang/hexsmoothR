# Get an appropriate UTM CRS for a study area

Determines the appropriate UTM coordinate reference system (CRS) for a
given study area based on the longitude/latitude of its centroid.

## Usage

``` r
get_utm_crs(study_area)
```

## Arguments

- study_area:

  sf object representing the study area. May be in any CRS; it will be
  transformed to WGS84 (EPSG:4326) internally for the centroid
  calculation.

## Value

Character string with the UTM CRS, e.g. \`"EPSG:32630"\` for UTM 30N.

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
study_area <- st_sf(geometry = st_sfc(
  st_polygon(list(matrix(c(-5, 35, 5, 35, 5, 45, -5, 45, -5, 35),
                         ncol = 2, byrow = TRUE))),
  crs = 4326
))
get_utm_crs(study_area)  # "EPSG:32630"
} # }
```
