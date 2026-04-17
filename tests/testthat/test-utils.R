test_that("ordinal_suffix handles teen and tens edge cases", {
  expect_equal(hexsmoothR:::ordinal_suffix(1),  "st")
  expect_equal(hexsmoothR:::ordinal_suffix(2),  "nd")
  expect_equal(hexsmoothR:::ordinal_suffix(3),  "rd")
  expect_equal(hexsmoothR:::ordinal_suffix(11), "th")
  expect_equal(hexsmoothR:::ordinal_suffix(12), "th")
  expect_equal(hexsmoothR:::ordinal_suffix(13), "th")
  expect_equal(hexsmoothR:::ordinal_suffix(21), "st")
  expect_equal(hexsmoothR:::ordinal_suffix(102), "nd")
})

test_that("hexsmoothR.verbose option suppresses console output", {
  withr::local_options(hexsmoothR.verbose = FALSE)
  out <- capture.output(hexsmoothR:::hex_msg("hello"))
  expect_length(out, 0L)
})

test_that("hex measurement helpers round-trip", {
  expect_equal(hex_edge_to_flat(hex_flat_to_edge(1234)), 1234)
  expect_equal(hex_circumradius_to_flat(hex_flat_to_circumradius(1234)), 1234)
})
