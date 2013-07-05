source("helper-tree.R")

context("Spline")

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
s <- new(Spline, FALSE)
test_that("Empty cubic splines behave sensibly", {
  expect_that(s$min, equals(Inf))
  expect_that(s$max, equals(-Inf))
  expect_that(s$eval(1), throws_error());
})

test_that("Spline contains correct data", {
  s$init(xx, yy)
  expect_that(s$xy,
              is_identical_to(cbind(xx, yy, deparse.level=0)))
  expect_that(c(s$min, s$max), is_identical_to(range(xx)))
})

test_that("Splines are accurate enough", {
  yy.C <- s$eval(xx.cmp)
  expect_that(yy.C, equals(yy.cmp, tolerance=1e-6))
})

## Akima splines:
a <- new(Spline, TRUE)
test_that("Empty Akima splines behave sensibly", {
  expect_that(a$min, equals(Inf))
  expect_that(a$max, equals(-Inf))
  expect_that(a$eval(1), throws_error());
})

test_that("Akima spline contains correct data", {
  a$init(xx, yy)
  expect_that(a$xy,
              is_identical_to(cbind(xx, yy, deparse.level=0)))
  expect_that(c(a$min, a$max), is_identical_to(range(xx)))
})

## Note that for the sin function, the Akima splines are way less
## accurate!
test_that("Akima splines are accurate enough", {
  yy.C <- a$eval(xx.cmp)
  expect_that(yy.C, equals(yy.cmp, tolerance=5e-5))
})
