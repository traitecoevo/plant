context("utils")

test_that("clamp_range", {
  f <- clamp_domain(identity, c(0, 1), Inf)
  expect_that(f(0), is_identical_to(0))
  expect_that(f(0L), is_identical_to(0.0))
  expect_that(f(-1), is_identical_to(Inf))
  expect_that(f(1.1), is_identical_to(Inf))

  f <- clamp_domain(identity, c(-Inf, 0), NA)
  expect_that(f(-100), equals(-100))

  expect_that(clamp_domain(identity, 1),
              throws_error("Expected length two range"))
  expect_that(clamp_domain(identity, c(0, 1), c(1, 1, 1)),
              throws_error("value must be length 1"))
})
