test_that("get_utm_crs auto-projects non-longlat input", {
  wgs <- test_polygon_wgs(-5, 35, 5, 45)
  utm_from_wgs <- get_utm_crs(wgs)
  utm_from_utm <- get_utm_crs(sf::st_transform(wgs, utm_from_wgs))
  expect_identical(utm_from_wgs, utm_from_utm)
})
