source("helper-tree.R")

context("Adaptive Quadrature (QAG)")

## First, on a simple function:
f <- function(x) sin(x)
g <- new(RFunctionWrapper, f, new.env())

a <- 0
b <- 1

ans.R <- integrate(f, a, b)
area.R <- ans.R$value
error.R <- ans.R$abs.error

test_that("R's integration behaves as expected", {
  expect_that(ans.R$subdivisions, equals(1))
})

test_that("Integration agrees with R on simple problem", {
  max.iterations <- 100
  eps <- .Machine$double.eps^0.25
  int.15 <- new(QAG, 15, max.iterations, eps, eps)
  expect_that(int.15$integrate(g, a, b), equals(area.R))
  expect_that(int.15$last_error,         equals(error.R))
  expect_that(int.15$last_iterations,    equals(1))
  int.21 <- new(QAG, 21, max.iterations, eps, eps)
  expect_that(int.21$integrate(g, a, b), equals(area.R))
  expect_that(int.21$last_error,         equals(error.R))
  expect_that(int.21$last_iterations,    equals(1))
  int.31 <- new(QAG, 31, max.iterations, eps, eps)
  expect_that(int.31$integrate(g, a, b), equals(area.R))
  expect_that(int.31$last_error,         equals(error.R))
  expect_that(int.31$last_iterations,    equals(1))
  int.41 <- new(QAG, 41, max.iterations, eps, eps)
  expect_that(int.41$integrate(g, a, b), equals(area.R))
  expect_that(int.41$last_error,         equals(error.R))
  expect_that(int.41$last_iterations,    equals(1))
  int.51 <- new(QAG, 51, max.iterations, eps, eps)
  expect_that(int.51$integrate(g, a, b), equals(area.R))
  expect_that(int.51$last_error,         equals(error.R))
  expect_that(int.51$last_iterations,    equals(1))
  int.61 <- new(QAG, 61, max.iterations, eps, eps)
  expect_that(int.61$integrate(g, a, b), equals(area.R))
  expect_that(int.61$last_error,         equals(error.R))
  expect_that(int.61$last_iterations,    equals(1))

  ## These should all be different.
  expect_that(identical(int.15$last_area, int.21$last_area),
              is_false())
  expect_that(identical(int.21$last_area, int.31$last_area),
              is_false())
  expect_that(identical(int.31$last_area, int.41$last_area),
              is_false())
  expect_that(identical(int.41$last_area, int.51$last_area),
              is_false())
  expect_that(identical(int.51$last_area, int.61$last_area),
              is_false())
})

## On a more complicated problem where subdivisions are required:

f <- function(x) 2^sin(sqrt(x))
g <- new(RFunctionWrapper, f, new.env())
a <- 0
b <- 50

ans.R <- integrate(f, a, b)
area.R <- ans.R$value
error.R <- ans.R$abs.error

test_that("R's integration behaves as expected", {
  expect_that(ans.R$subdivisions, is_greater_than(1))
})

test_that("Integration agrees with R on complex problem", {
  ## It looks like the error is computed slightly differently between
  ## these...
  max.iterations <- 100
  eps <- .Machine$double.eps^0.25
  int.15 <- new(QAG, 15, max.iterations, eps, eps)
  expect_that(int.15$integrate(g, a, b), equals(area.R, tolerance=eps))
  ## expect_that(int.15$last_error,         equals(error.R))
  expect_that(int.15$last_iterations,    is_greater_than(1))
  int.21 <- new(QAG, 21, max.iterations, eps, eps)
  expect_that(int.21$integrate(g, a, b), equals(area.R, tolerance=eps))
  ## expect_that(int.21$last_error,         equals(error.R))
  expect_that(int.21$last_iterations,    is_greater_than(1))
  int.31 <- new(QAG, 31, max.iterations, eps, eps)
  expect_that(int.31$integrate(g, a, b), equals(area.R, tolerance=eps))
  ## expect_that(int.31$last_error,         equals(error.R))
  expect_that(int.31$last_iterations,    is_greater_than(1))
  int.41 <- new(QAG, 41, max.iterations, eps, eps)
  expect_that(int.41$integrate(g, a, b), equals(area.R, tolerance=eps))
  ## expect_that(int.41$last_error,         equals(error.R))
  expect_that(int.41$last_iterations,    is_greater_than(1))
  int.51 <- new(QAG, 51, max.iterations, eps, eps)
  expect_that(int.51$integrate(g, a, b), equals(area.R, tolerance=eps))
  ## expect_that(int.51$last_error,         equals(error.R))
  expect_that(int.51$last_iterations,    is_greater_than(1))
  int.61 <- new(QAG, 61, max.iterations, eps, eps)
  expect_that(int.61$integrate(g, a, b), equals(area.R, tolerance=eps))
  ## expect_that(int.61$last_error,         equals(error.R))
  expect_that(int.61$last_iterations,    is_greater_than(1))

  ## These should all be different.
  expect_that(identical(int.15$last_area, int.21$last_area),
              is_false())
  expect_that(identical(int.21$last_area, int.31$last_area),
              is_false())
  expect_that(identical(int.31$last_area, int.41$last_area),
              is_false())
  expect_that(identical(int.41$last_area, int.51$last_area),
              is_false())
  expect_that(identical(int.51$last_area, int.61$last_area),
              is_false())
})

test_that("Cannot make non-existent rules", {
  max.iterations <- 100
  eps <- 1e-8
  expect_that(new(QAG, 0,   max.iterations, eps, eps), throws_error())
  expect_that(new(QAG, 100, max.iterations, eps, eps), throws_error())
  expect_that(new(QAG, NA,  max.iterations, eps, eps), throws_error())
})

test_that("Remembering intervals works", {
  max.iterations <- 100
  eps <- .Machine$double.eps^0.25

  int.15 <- new(QAG, 15, max.iterations, eps, eps)
  ans.15 <- int.15$integrate(g, a, b)
  err.15 <- int.15$last_error
  expect_that(ans.15, equals(area.R, tolerance=eps))

  intervals <- int.15$last_intervals
  tmp <- do.call(cbind, intervals)
  tmp <- tmp[order(tmp[,1]),]

  ## check that the intervals make sense:
  expect_that(tmp[1,1], is_identical_to(a))
  expect_that(tmp[nrow(tmp),2], is_identical_to(b))
  expect_that(tmp[-1,1], is_identical_to(tmp[-nrow(tmp),2]))
  expect_that(all(diff(tmp[,1]) > 0), is_true())

  ## Now, do the integration with a new set of intervals:
  mid <- (a + b)/2
  intervals.2 <- list(c(a, mid), c(mid, b))
  ans.2 <- int.15$integrate_with_intervals(g, intervals.2)
  err.2 <- int.15$last_error
  expect_that(err.2, is_greater_than(err.15))
  expect_that(identical(err.15, err.2), is_false())

  ## Then with the computed intervals;
  ans.redo <- int.15$integrate_with_intervals(g, intervals)
  err.redo <- int.15$last_error
  expect_that(ans.redo, is_identical_to(ans.15))
  expect_that(err.redo, equals(err.15))
})

test_that("Non-adaptive integration works", {
  max.iterations <- 100
  eps <- .Machine$double.eps^0.25

  for (rule in c(15, 61)) {
    int.a <- new(QAG, rule, max.iterations, eps, eps)
    int.f <- new(QAG, rule)
    int.q <- new(QK,  rule)

    ans.a <- int.a$integrate(g, a, b)
    ans.f <- int.f$integrate(g, a, b)
    ans.q <- int.q$integrate(g, a, b)

    expect_that(ans.f, is_identical_to(ans.q))
    expect_that(int.f$last_error,
                is_identical_to(int.q$last_error))
  }
})
