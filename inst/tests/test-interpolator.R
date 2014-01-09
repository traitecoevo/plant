source("helper-tree.R")

context("Interpolator")

target <- function(x) sin(2*x)

xx <- seq(0, 2*pi, length=101)
yy <- target(xx)

## Places to evaluate, and the true values at that point:
xx.cmp <- seq(min(xx), max(xx), length=2*length(xx))
yy.cmp <- target(xx.cmp)

test_that("Test data is sensible", {
  yy.R <- spline(xx, yy, xout=xx.cmp)$y
  expect_that(yy.R, equals(yy.cmp, tolerance=1e-6))
})

## Basic cubic splines:
s <- new(Interpolator, FALSE, FALSE)
test_that("Empty cubic splines behave sensibly", {
  expect_that(s$min, equals(Inf))
  expect_that(s$max, equals(-Inf))
  expect_that(s$eval(1), throws_error());
  expect_that(s$type,    throws_error()); #  stupid
})

test_that("Splines require sensible data", {
  expect_that(s$init(c(1, 2, 3), 1), throws_error())
  expect_that(s$init(numeric(0), numeric(0)), throws_error())
})

test_that("Spline contains correct data", {
  s$init(xx, yy)
  expect_that(s$xy,
              is_identical_to(cbind(xx, yy, deparse.level=0)))
  expect_that(c(s$min, s$max), is_identical_to(range(xx)))
  expect_that(s$type, equals("cspline"))
})

test_that("Splines are accurate enough", {
  yy.C <- s$eval(xx.cmp)
  expect_that(yy.C, equals(yy.cmp, tolerance=1e-6))
})

test_that("Spline derivatives are correct", {
  cmp.deriv <- splinefun(xx, yy)(xx.cmp, deriv=1L)
  expect_that(s$deriv(xx.cmp), equals(cmp.deriv, tolerance=1e-5))
})

## Akima splines:
a <- new(Interpolator, TRUE, FALSE)
test_that("Empty Akima splines behave sensibly", {
  expect_that(a$min, equals(Inf))
  expect_that(a$max, equals(-Inf))
  expect_that(a$eval(1), throws_error());
  expect_that(a$type,    throws_error()); #  stupid
})

test_that("Akima spline contains correct data", {
  a$init(xx, yy)
  expect_that(a$xy,
              is_identical_to(cbind(xx, yy, deparse.level=0)))
  expect_that(c(a$min, a$max), is_identical_to(range(xx)))
  expect_that(a$type,  equals("akima"))
})

## Note that for the sin function, the Akima splines are way less
## accurate!
test_that("Akima splines are accurate enough", {
  yy.C <- a$eval(xx.cmp)
  expect_that(yy.C, equals(yy.cmp, tolerance=5e-5))
})

## Linear interpolation
l <- new(Interpolator, FALSE, TRUE)
test_that("Empty linear interpolation behave sensibly", {
  expect_that(l$min, equals(Inf))
  expect_that(l$max, equals(-Inf))
  expect_that(l$eval(1), throws_error());
  expect_that(l$type,    throws_error()); #  stupid
})

test_that("Linear interpolation contains correct data", {
  l$init(xx, yy)
  expect_that(l$xy,
              is_identical_to(cbind(xx, yy, deparse.level=0)))
  expect_that(c(l$min, l$max), is_identical_to(range(xx)))
  expect_that(l$type, equals("linear"))
})

test_that("Linear interpolation is accurate enough", {
  yy.C <- l$eval(xx.cmp)
  yy.cmp <- approx(xx, yy, xout=xx.cmp)$y
  expect_that(yy.C, equals(yy.cmp))
})
