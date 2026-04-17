test_that("BFS neighbours: order 1 and order 2 are disjoint and non-empty", {
  grid <- test_hex_grid()
  topo <- suppressMessages(compute_topology(grid, neighbor_orders = 2))

  o1 <- topo$neighbors[[1]]
  o2 <- topo$neighbors[[2]]

  expect_true(length(o1) == nrow(grid))
  expect_true(length(o2) == nrow(grid))

  for (i in seq_along(o1)) {
    expect_length(intersect(o1[[i]], o2[[i]]), 0L)
    expect_false(i %in% o1[[i]])
    expect_false(i %in% o2[[i]])
  }
})

test_that("a list of grids returns a named list of topologies", {
  grid <- test_hex_grid()
  out <- suppressMessages(
    compute_topology(list(a = grid, b = grid), neighbor_orders = 2)
  )

  expect_named(out, c("a", "b"))
  expect_true(all(vapply(out, is.list, logical(1))))
})
