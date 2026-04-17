test_that("'raw' equals the original (unsmoothed) value", {
  grid <- test_hex_grid()
  topo <- suppressMessages(compute_topology(grid))
  vals <- runif(nrow(grid), 1, 10)

  out <- suppressMessages(suppressWarnings(
    smooth_variables(
      variable_values = list(v = vals),
      neighbors = topo$neighbors,
      weights   = topo$weights,
      var_names = "v"
    )
  ))

  expect_equal(out$v$raw, vals)
  expect_false(isTRUE(all.equal(out$v$raw, out$v$weighted_combined)))
})

test_that("input lengths are validated", {
  grid <- test_hex_grid()
  topo <- suppressMessages(compute_topology(grid))

  expect_error(
    suppressMessages(smooth_variables(
      variable_values = list(a = 1:3, b = 1:4),
      neighbors = topo$neighbors,
      weights   = topo$weights,
      var_names = c("a", "b")
    )),
    regexp = "same length"
  )
})

test_that("var_names are inferred from list names", {
  grid <- test_hex_grid()
  topo <- suppressMessages(compute_topology(grid))
  vals <- list(foo = runif(nrow(grid), 1, 10))

  out <- suppressMessages(suppressWarnings(
    smooth_variables(
      variable_values = vals,
      neighbors = topo$neighbors,
      weights   = topo$weights
    )
  ))

  expect_named(out, "foo")
})
