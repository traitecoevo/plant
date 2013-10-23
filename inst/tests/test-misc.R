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

## Integration works?

test_that("GSL integration rule translation correct", {
  ## From string -> key
  expect_that(integrator_gsl_rule("GAUSS15"), equals(1))
  expect_that(integrator_gsl_rule("GAUSS21"), equals(2))
  expect_that(integrator_gsl_rule("GAUSS31"), equals(3))
  expect_that(integrator_gsl_rule("GAUSS41"), equals(4))
  expect_that(integrator_gsl_rule("GAUSS51"), equals(5))
  expect_that(integrator_gsl_rule("GAUSS61"), equals(6))
  expect_that(integrator_gsl_rule("QAGS"),        equals(-1))
  expect_that(integrator_gsl_rule("NONADAPTIVE"), equals(-2))
  expect_that(integrator_gsl_rule("nonexistant"), throws_error())

  ## and from key -> string
  expect_that(integrator_gsl_rule_name( 1), equals("GAUSS15"))
  expect_that(integrator_gsl_rule_name( 2), equals("GAUSS21"))
  expect_that(integrator_gsl_rule_name( 3), equals("GAUSS31"))
  expect_that(integrator_gsl_rule_name( 4), equals("GAUSS41"))
  expect_that(integrator_gsl_rule_name( 5), equals("GAUSS51"))
  expect_that(integrator_gsl_rule_name( 6), equals("GAUSS61"))
  expect_that(integrator_gsl_rule_name(-1), equals("QAGS"))
  expect_that(integrator_gsl_rule_name(-2), equals("NONADAPTIVE"))
  expect_that(integrator_gsl_rule_name( 7), throws_error())
})

test_that("Integration works", {
  f <- function(rule)
    test_integrator(pars, x.min, x.max, rule)
  f.int <- diff(with(as.list(pars), function(x)
                     a/3*x^3 + b/2*x^2 + c*x)(c(x.min, x.max)))

  expect_that(f("GAUSS15"), equals(f.int))
  expect_that(f("GAUSS21"), equals(f.int))
  expect_that(f("GAUSS31"), equals(f.int))
  expect_that(f("GAUSS41"), equals(f.int))
  expect_that(f("GAUSS51"), equals(f.int))
  expect_that(f("GAUSS61"), equals(f.int))
  expect_that(f("QAGS"),        equals(f.int))
  expect_that(f("NONADAPTIVE"), equals(f.int))

  expect_that(f("nonexistant"), throws_error())
  ## Different algorithms differ in answer slightly (not all do!)
  expect_that(f("QAGS")    != f("GAUSS15"), is_true())
  expect_that(f("GAUSS51") != f("GAUSS15"), is_true())
})

set.seed(1)
a <- runif(10)
b <- runif(length(a))
expect_that(test_sum_double(a, b), is_identical_to(a + b))

a <- as.integer(a * 10)
b <- as.integer(b * 10)
expect_that(test_sum_int(a, b), is_identical_to(a + b))

expect_that(test_to_rcpp_matrix(list(a, b)),
            is_identical_to(cbind(a, b, deparse.level=0)))
expect_that(test_from_rcpp_matrix(cbind(a, b)),
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
