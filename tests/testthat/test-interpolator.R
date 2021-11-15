context("Interpolator")

target <- function(x) sin(2*x)

xx <- seq(0, 2*pi, length.out=101)
yy <- target(xx)

## Places to evaluate, and the true values at that point:
xx_cmp <- seq(min(xx), max(xx), length.out=2*length(xx))
yy_cmp <- target(xx_cmp)

test_that("Test data is sensible", {
  yy_R <- spline(xx, yy, xout=xx_cmp)$y
  expect_equal(yy_R, yy_cmp, tolerance=1e-6)
})

## Basic cubic splines:
test_that("Empty cubic splines behave sensibly", {
  s <- Interpolator()
  expect_equal(s$min, Inf)
  expect_equal(s$max, -Inf)
  expect_error(s$eval(1), "Interpolator not initialised")
})

test_that("Splines require sensible data", {
  s <- Interpolator()
  expect_error(s$init(c(1, 2, 3), 1), "Incorrect length input")
  expect_error(s$init(numeric(0), numeric(0)), "insufficient number of points")
  expect_error(s$init(c(3, 2, 1), c(1, 1, 1)), "spline control points must be in non-descending order")
})

test_that("Spline contains correct data", {
  s <- Interpolator()
  s$init(xx, yy)
  expect_identical(s$xy, cbind(xx, yy, deparse.level=0))
  expect_identical(c(s$min, s$max), range(xx))
})

test_that("Splines are accurate enough", {
  s <- Interpolator()
  s$init(xx, yy)
  yy_C <- s$eval(xx_cmp)
  expect_equal(yy_C, yy_cmp, tolerance=1e-6)
})
