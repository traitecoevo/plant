if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("Interpolator")

target <- function(x) sin(2*x)

xx <- seq(0, 2*pi, length=101)
yy <- target(xx)

## Places to evaluate, and the true values at that point:
xx_cmp <- seq(min(xx), max(xx), length=2*length(xx))
yy_cmp <- target(xx_cmp)

test_that("Test data is sensible", {
  yy_R <- spline(xx, yy, xout=xx_cmp)$y
  expect_that(yy_R, equals(yy_cmp, tolerance=1e-6))
})

## Basic cubic splines:
test_that("Empty cubic splines behave sensibly", {
  s <- Interpolator()
  expect_that(s$min, equals(Inf))
  expect_that(s$max, equals(-Inf))
  expect_that(s$eval(1),
              throws_error("Interpolator not initialised"))
  expect_that(s$type,
              throws_error("Interpolator not initialised")) #  stupid
})

test_that("Splines require sensible data", {
  s <- Interpolator()
  expect_that(s$init(c(1, 2, 3), 1),
              throws_error("Incorrect length input"))
  expect_that(s$init(numeric(0), numeric(0)),
              throws_error("insufficient number of points"))
})

test_that("Spline contains correct data", {
  s <- Interpolator()
  s$init(xx, yy)
  expect_that(s$xy,
              is_identical_to(cbind(xx, yy, deparse.level=0)))
  expect_that(c(s$min, s$max), is_identical_to(range(xx)))
  expect_that(s$type, equals("cspline"))
})

test_that("Splines are accurate enough", {
  s <- Interpolator()
  s$init(xx, yy)
  yy_C <- s$eval(xx_cmp)
  expect_that(yy_C, equals(yy_cmp, tolerance=1e-6))

  s2 <- Interpolator_TK()
  s2$init(xx, yy)
  yy_C2 <- s2$eval(xx_cmp)
  expect_that(yy_C2, equals(yy_cmp, tolerance=1e-6))
  expect_that(yy_C2, not(is_identical_to(yy_C)))
})

test_that("Spline derivatives are correct", {
  s <- Interpolator()
  s$init(xx, yy)
  cmp_deriv <- splinefun(xx, yy)(xx_cmp, deriv=1L)
  expect_that(s$deriv(xx_cmp), equals(cmp_deriv, tolerance=1e-5))
})
