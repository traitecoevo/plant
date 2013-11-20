source("helper-tree.R")

context("Misc")

pars <- c(a=1, b=0, c=-5)
x.min <- 0
x.max <- 5
xx <- c(x.min, 2, x.max)
f <- with(as.list(pars), function(x) a*x*x + b*x + c)

## Functor evaluates correctly?
expect_that(test_functor(xx, pars),
            equals(f(xx)))

## Test with RFunctionWrapper:
f.wrapped <- new(RFunctionWrapper, f, new.env())
expect_that(sapply(xx, function(x) f.wrapped$target(x)),
            is_identical_to(f(xx)))

## Find root correctly?
## uniroot(f, c(x.min, x.max))$root

## Analytically
tmp <- with(as.list(pars), (-b + c(-1,1) * sqrt(b^2 - 4*a*c)/(2*a)))
f.root <- tmp[tmp > x.min & tmp < x.max]

## quick sanity check.
expect_that(length(f.root), equals(1))

expect_that(test_find_root(pars, x.min, x.max),
            equals(f.root, tolerance=1e-5))

## Root finding to particular value
value <- 10.0
f.value <- uniroot(function(x) f(x) - value, c(x.min, x.max))$root
expect_that(test_find_value(pars, value, x.min, x.max),
            equals(f.value, tolerance=1e-5))

## What happens if we're out of range?
expect_that(test_find_value(pars, value, x.min, x.min+1),
            throws_error())

set.seed(1)
a <- runif(10)
b <- runif(length(a))
expect_that(test_sum_double(a, b), is_identical_to(a + b))

expect_that(test_to_rcpp_numeric_matrix(list(a, b)),
            is_identical_to(cbind(a, b, deparse.level=0)))
expect_that(test_from_rcpp_numeric_matrix(cbind(a, b)),
            is_identical_to(list(a, b)))

a <- as.integer(a * 10)
b <- as.integer(b * 10)
expect_that(test_sum_int(a, b), is_identical_to(a + b))

expect_that(test_to_rcpp_integer_matrix(list(a, b)),
            is_identical_to(cbind(a, b, deparse.level=0)))
expect_that(test_from_rcpp_integer_matrix(cbind(a, b)),
            is_identical_to(list(a, b)))

## Test the simple finite differencing gradient function.
gradient.fd.forward <- function(f, x, dx)
  (f(x + dx) - f(x)) / dx
gradient.fd.centre <- function(f, x, dx)
  (f(x + dx/2) - f(x - dx/2)) / dx
gradient.fd.backward <- function(f, x, dx)
  (f(x - dx) - f(x)) / (-dx)

dx <- 0.001
x <- 1
expect_that(test_gradient(x, dx, 1, pars),
            is_identical_to(gradient.fd.forward(f, x, dx)))
expect_that(test_gradient(x, dx, 0, pars),
            is_identical_to(gradient.fd.centre(f, x, dx)))
expect_that(test_gradient(x, dx, -1, pars),
            is_identical_to(gradient.fd.backward(f, x, dx)))

library(numDeriv)
method.args <- list(d=1e-6, eps=1e-6)
expect_that(test_gradient_richardson(x, 1e-6, 4L, pars),
            is_identical_to(grad(f, x, method.args=method.args)))

n <- 20
set.seed(1)
xx <- sort(runif(n))
yy <- runif(n)
expect_that(trapezium(xx, yy),
            equals(sum((xx[-1] - xx[-n]) * (yy[-1] + yy[-n])) / 2))
expect_that(trapezium(c(1, 1), c(1, 2)),
            is_identical_to(0.0))

test_that("Trapezium local error estimate is correct", {
  ## Here is the naive implementation of integration error...
  local.error.integration <- function(x, y) {
    n <- length(x)
    i0 <- -c(n-1, n)
    i1 <- -c(1,   n)
    i2 <- -c(1,   2)

    y1.fine <-
      0.5*((y[i0] + y[i1])*(x[i1] - x[i0]) +
           (y[i1] + y[i2])*(x[i2] - x[i1]))
    y1.coarse <-
      0.5*((y[i0] + y[i2])*(x[i2] - x[i0]))
    tot <- abs(trapezium(x, y))
    y1.fine   <- y1.fine / tot
    y1.coarse <- y1.coarse / tot

    err.abs <- abs(    y1.coarse - y1.fine)
    err.rel <- abs(1 - y1.coarse / y1.fine)
    err.rel[is.nan(err.rel)] <- Inf
    c(NA, pmin(err.abs, err.rel), NA)
  }

  ## Some data:
  set.seed(1)
  xx <- sort(runif(100, 0, 6*pi))
  yy <- sin(xx) + 1

  expect_that(local_error_integration(xx, yy),
              equals(local.error.integration(xx, yy)))

  expect_that(local_error_integration(numeric(0), numeric(0)),
              equals(numeric(0)))
  expect_that(local_error_integration(xx[1], yy[1]),
              equals(NA_real_))
  expect_that(local_error_integration(xx[1:2], yy[1:2]),
              equals(c(NA_real_, NA_real_)))
  expect_that(local_error_integration(xx[1:2], yy[1:3]),
              throws_error())
})
