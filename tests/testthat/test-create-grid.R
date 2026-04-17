test_that("create_grid auto-detects projection_crs when NULL", {
  wgs <- test_polygon_wgs(-5, 35, 5, 45)
  expect_no_error(
    suppressMessages(create_grid(wgs, cell_size = 20000, type = "hexagonal"))
  )
})
