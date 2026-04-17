# Shared test fixtures.
# `helper-*.R` files are sourced automatically by testthat before each test file.

suppressPackageStartupMessages({
  library(sf)
})

#' Build a small WGS84 polygon for tests.
#' @keywords internal
test_polygon_wgs <- function(xmin = 0, ymin = 0, xmax = 1, ymax = 1) {
  sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(matrix(
      c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
      ncol = 2, byrow = TRUE
    ))),
    crs = 4326
  ))
}

#' Build a small UTM hexagonal grid suitable for fast smoothing tests.
#' @keywords internal
test_hex_grid <- function(cell_size = 20000) {
  wgs <- test_polygon_wgs()
  utm <- get_utm_crs(wgs)
  utm_area <- sf::st_transform(wgs, utm)
  suppressMessages(
    create_grid(utm_area, cell_size = cell_size, type = "hexagonal")
  )
}
