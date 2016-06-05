context("utils")

test_that("clamp_range", {
  f <- clamp_domain(identity, c(0, 1), Inf)
  expect_identical(f(0), 0)
  expect_identical(f(0L), 0.0)
  expect_identical(f(-1), Inf)
  expect_identical(f(1.1), Inf)

  f <- clamp_domain(identity, c(-Inf, 0), NA)
  expect_equal(f(-100), -100)

  expect_error(clamp_domain(identity, 1), "Expected length two range")
  expect_error(clamp_domain(identity, c(0, 1), c(1, 1, 1)), "value must be length 1")
})
