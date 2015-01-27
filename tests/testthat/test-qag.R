if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("QAG")

## First where no subdivisions are required:

test_that("Integration agrees with R on simple problem", {
  f <- sin
  a <- 0
  b <- 1

  ans_R <- integrate(f, a, b)
  expect_that(ans_R$subdivisions, equals(1))

  area_R <- ans_R$value
  error_R <- ans_R$abs.error

  max_iterations <- 100
  eps <- .Machine$double.eps^0.25

  int_15 <- QAG(15, max_iterations, eps, eps)
  expect_that(int_15$integrate(f, a, b), equals(area_R))
  expect_that(int_15$last_error,         equals(error_R))
  expect_that(int_15$last_iterations,    equals(1))
  int_21 <- QAG(21, max_iterations, eps, eps)
  expect_that(int_21$integrate(f, a, b), equals(area_R))
  expect_that(int_21$last_error,         equals(error_R))
  expect_that(int_21$last_iterations,    equals(1))
  int_31 <- QAG(31, max_iterations, eps, eps)
  expect_that(int_31$integrate(f, a, b), equals(area_R))
  expect_that(int_31$last_error,         equals(error_R))
  expect_that(int_31$last_iterations,    equals(1))
  int_41 <- QAG(41, max_iterations, eps, eps)
  expect_that(int_41$integrate(f, a, b), equals(area_R))
  expect_that(int_41$last_error,         equals(error_R))
  expect_that(int_41$last_iterations,    equals(1))
  int_51 <- QAG(51, max_iterations, eps, eps)
  expect_that(int_51$integrate(f, a, b), equals(area_R))
  expect_that(int_51$last_error,         equals(error_R))
  expect_that(int_51$last_iterations,    equals(1))
  int_61 <- QAG(61, max_iterations, eps, eps)
  expect_that(int_61$integrate(f, a, b), equals(area_R))
  expect_that(int_61$last_error,         equals(error_R))
  expect_that(int_61$last_iterations,    equals(1))

  ## These should all be different, though the first case is not
  ## necessarily so.
  ## expect_that(identical(int_15$last_area, int_21$last_area),
  ##            is_false())
  expect_that(identical(int_21$last_area, int_31$last_area),
              is_false())
  expect_that(identical(int_31$last_area, int_41$last_area),
              is_false())
  expect_that(identical(int_41$last_area, int_51$last_area),
              is_false())
  expect_that(identical(int_51$last_area, int_61$last_area),
              is_false())
})

test_that("Integration agrees when subdividing", {
  f <- function(x) 2^sin(sqrt(x))
  a <- 0
  b <- 50

  ans_R <- integrate(f, a, b)
  area_R <- ans_R$value
  error_R <- ans_R$abs.error

  expect_that(ans_R$subdivisions, is_greater_than(1))

  max_iterations <- 100
  eps <- .Machine$double.eps^0.25

  int_15 <- QAG(15, max_iterations, eps, eps)
  expect_that(int_15$integrate(f, a, b), equals(area_R, tolerance=eps))
  ## expect_that(int_15$last_error,         equals(error_R))
  expect_that(int_15$last_iterations,    is_greater_than(1))
  int_21 <- QAG(21, max_iterations, eps, eps)
  expect_that(int_21$integrate(f, a, b), equals(area_R, tolerance=eps))
  ## expect_that(int_21$last_error,         equals(error_R))
  expect_that(int_21$last_iterations,    is_greater_than(1))
  int_31 <- QAG(31, max_iterations, eps, eps)
  expect_that(int_31$integrate(f, a, b), equals(area_R, tolerance=eps))
  ## expect_that(int_31$last_error,         equals(error_R))
  expect_that(int_31$last_iterations,    is_greater_than(1))
  int_41 <- QAG(41, max_iterations, eps, eps)
  expect_that(int_41$integrate(f, a, b), equals(area_R, tolerance=eps))
  ## expect_that(int_41$last_error,         equals(error_R))
  expect_that(int_41$last_iterations,    is_greater_than(1))
  int_51 <- QAG(51, max_iterations, eps, eps)
  expect_that(int_51$integrate(f, a, b), equals(area_R, tolerance=eps))
  ## expect_that(int_51$last_error,         equals(error_R))
  expect_that(int_51$last_iterations,    is_greater_than(1))
  int_61 <- QAG(61, max_iterations, eps, eps)
  expect_that(int_61$integrate(f, a, b), equals(area_R, tolerance=eps))
  ## expect_that(int_61$last_error,         equals(error_R))
  expect_that(int_61$last_iterations,    is_greater_than(1))

  ## These should all be different.
  expect_that(identical(int_15$last_area, int_21$last_area),
              is_false())
  expect_that(identical(int_21$last_area, int_31$last_area),
              is_false())
  expect_that(identical(int_31$last_area, int_41$last_area),
              is_false())
  expect_that(identical(int_41$last_area, int_51$last_area),
              is_false())
  expect_that(identical(int_51$last_area, int_61$last_area),
              is_false())
})

test_that("Cannot make non-existent rules", {
  max_iterations <- 100
  eps <- 1e-8
  expect_that(QAG(0,   max_iterations, eps, eps),
              throws_error("Unknown rule 0"))
  expect_that(QAG(100, max_iterations, eps, eps),
              throws_error("Unknown rule 100"))
  expect_that(QAG(NA,  max_iterations, eps, eps),
              throws_error("Unknown rule"))
})

test_that("Remembering intervals works", {
  f <- function(x) 2^sin(sqrt(x))
  a <- 0
  b <- 50
  max_iterations <- 100
  eps <- .Machine$double.eps^0.25

  ans_R <- integrate(f, a, b)
  area_R <- ans_R$value
  error_R <- ans_R$abs.error

  int_15 <- QAG(15, max_iterations, eps, eps)
  ## Can't reuse intervals before they're stored:
  expect_that(int_15$integrate_with_last_intervals(f, a, b),
              throws_error("No stored intervals to use"))
  ans_15 <- int_15$integrate(f, a, b)
  err_15 <- int_15$last_error
  expect_that(ans_15, equals(area_R, tolerance=eps))

  intervals <- int_15$last_intervals
  tmp <- do.call(cbind, intervals)
  tmp <- tmp[order(tmp[,1]),]

  ## check that the intervals make sense:
  expect_that(tmp[1,1], is_identical_to(a))
  expect_that(tmp[nrow(tmp),2], is_identical_to(b))
  expect_that(tmp[-1,1], is_identical_to(tmp[-nrow(tmp),2]))
  expect_that(all(diff(tmp[,1]) > 0), is_true())

  ## Reusing the last intervals is easy:
  expect_that(int_15$integrate_with_last_intervals(f, a, b),
              is_identical_to(ans_15))
  expect_that(int_15$last_intervals, is_identical_to(intervals))

  ## Now, do the integration with a new set of intervals:
  mid <- (a + b)/2
  intervals_2 <- list(c(a, mid), c(mid, b))
  ans_2 <- int_15$integrate_with_intervals(f, intervals_2)
  err_2 <- int_15$last_error
  expect_that(err_2, is_greater_than(err_15))
  expect_that(identical(err_15, err_2), is_false())
  expect_that(int_15$last_intervals,
              equals(intervals_2))

  ## Then with the computed intervals;
  ans_redo <- int_15$integrate_with_intervals(f, intervals)
  err_redo <- int_15$last_error
  expect_that(ans_redo, is_identical_to(ans_15))
  expect_that(err_redo, equals(err_15))
  expect_that(int_15$last_intervals, is_identical_to(intervals))

  ## Rescaling intervals is hard to test well, but here goes:
  b1 <- 49.5
  ans_scal <- int_15$integrate_with_last_intervals(f, 0, b1)

  intervals_scal <- int_15$last_intervals
  tmp_scal <- do.call(cbind, intervals_scal)
  tmp_scal <- tmp_scal[order(tmp_scal[,1]),]
  expect_that(dim(tmp_scal), is_identical_to(dim(tmp)))
  expect_that(tmp_scal, equals(tmp * b1 / b))
})

test_that("Non-adaptive integration works", {
  f <- function(x) 2^sin(sqrt(x))
  a <- 0
  b <- 50
  max_iterations <- 100
  eps <- .Machine$double.eps^0.25

  for (rule in range(qk_rules())) {
    int_a <- QAG(rule, max_iterations, eps, eps)
    int_f <- QAG(rule, 1,              NA,  NA)
    int_q <- QK(rule)

    expect_that(int_a$is_adaptive, is_true())
    expect_that(int_f$is_adaptive, is_false())

    ans_a <- int_a$integrate(f, a, b)
    ans_f <- int_f$integrate(f, a, b)
    ans_q <- int_q$integrate(f, a, b)

    expect_that(ans_f, is_identical_to(ans_q))
    expect_that(int_f$last_error,
                is_identical_to(int_q$last_error))
    if (length(int_a$last_intervals[[1]]) > 1L) {
      expect_that(ans_a, not(is_identical_to(ans_q)))
    } else {
      expect_that(ans_a, is_identical_to(ans_q))
    }
  }
})
