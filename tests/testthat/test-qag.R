context("QAG")

## First where no subdivisions are required:

test_that("Integration agrees with R on simple problem", {
  f <- sin
  a <- 0
  b <- 1

  ans_R <- integrate(f, a, b)
  expect_equal(ans_R$subdivisions, 1)

  area_R <- ans_R$value
  error_R <- ans_R$abs.error

  max_iterations <- 100
  eps <- .Machine$double.eps^0.25

  int_15 <- QAG(15, max_iterations, eps, eps)
  expect_equal(int_15$integrate(f, a, b), area_R)
  expect_equal(int_15$last_error, error_R)
  expect_equal(int_15$last_iterations, 1)
  int_21 <- QAG(21, max_iterations, eps, eps)
  expect_equal(int_21$integrate(f, a, b), area_R)
  expect_equal(int_21$last_error, error_R)
  expect_equal(int_21$last_iterations, 1)
  int_31 <- QAG(31, max_iterations, eps, eps)
  expect_equal(int_31$integrate(f, a, b), area_R)
  expect_equal(int_31$last_error, error_R)
  expect_equal(int_31$last_iterations, 1)
  int_41 <- QAG(41, max_iterations, eps, eps)
  expect_equal(int_41$integrate(f, a, b), area_R)
  expect_equal(int_41$last_error, error_R)
  expect_equal(int_41$last_iterations, 1)
  int_51 <- QAG(51, max_iterations, eps, eps)
  expect_equal(int_51$integrate(f, a, b), area_R)
  expect_equal(int_51$last_error, error_R)
  expect_equal(int_51$last_iterations, 1)
  int_61 <- QAG(61, max_iterations, eps, eps)
  expect_equal(int_61$integrate(f, a, b), area_R)
  expect_equal(int_61$last_error, error_R)
  expect_equal(int_61$last_iterations, 1)

  ## These should all be different, though the first case is not
  ## necessarily so.
  ## expect_false(identical(int_15$last_area, int_21$last_area))
  expect_false(identical(int_21$last_area, int_31$last_area))
  expect_false(identical(int_31$last_area, int_41$last_area))
  expect_false(identical(int_41$last_area, int_51$last_area))
  expect_false(identical(int_51$last_area, int_61$last_area))
})

test_that("Integration agrees when subdividing", {
  f <- function(x) 2^sin(sqrt(x))
  a <- 0
  b <- 50

  ans_R <- integrate(f, a, b)
  area_R <- ans_R$value
  error_R <- ans_R$abs.error

  expect_gt(ans_R$subdivisions, 1)

  max_iterations <- 100
  eps <- .Machine$double.eps^0.25

  int_15 <- QAG(15, max_iterations, eps, eps)
  expect_equal(int_15$integrate(f, a, b), area_R, tolerance=eps)
  ## expect_equal(int_15$last_error, error_R)
  expect_gt(int_15$last_iterations, 1)
  int_21 <- QAG(21, max_iterations, eps, eps)
  expect_equal(int_21$integrate(f, a, b), area_R, tolerance=eps)
  ## expect_equal(int_21$last_error, error_R)
  expect_gt(int_21$last_iterations, 1)
  int_31 <- QAG(31, max_iterations, eps, eps)
  expect_equal(int_31$integrate(f, a, b), area_R, tolerance=eps)
  ## expect_equal(int_31$last_error, error_R)
  expect_gt(int_31$last_iterations, 1)
  int_41 <- QAG(41, max_iterations, eps, eps)
  expect_equal(int_41$integrate(f, a, b), area_R, tolerance=eps)
  ## expect_equal(int_41$last_error, error_R)
  expect_gt(int_41$last_iterations, 1)
  int_51 <- QAG(51, max_iterations, eps, eps)
  expect_equal(int_51$integrate(f, a, b), area_R, tolerance=eps)
  ## expect_equal(int_51$last_error, error_R)
  expect_gt(int_51$last_iterations, 1)
  int_61 <- QAG(61, max_iterations, eps, eps)
  expect_equal(int_61$integrate(f, a, b), area_R, tolerance=eps)
  ## expect_equal(int_61$last_error, error_R)
  expect_gt(int_61$last_iterations, 1)

  ## These should all be different.
  expect_false(identical(int_15$last_area, int_21$last_area))
  expect_false(identical(int_21$last_area, int_31$last_area))
  expect_false(identical(int_31$last_area, int_41$last_area))
  expect_false(identical(int_41$last_area, int_51$last_area))
  expect_false(identical(int_51$last_area, int_61$last_area))
})

test_that("Cannot make non-existent rules", {
  max_iterations <- 100
  eps <- 1e-8
  expect_error(QAG(0,   max_iterations, eps, eps), "Unknown rule 0")
  expect_error(QAG(100, max_iterations, eps, eps), "Unknown rule 100")
  expect_error(QAG(NA,  max_iterations, eps, eps), "Unknown rule")
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
  expect_error(int_15$integrate_with_last_intervals(f, a, b), "No stored intervals to use")
  ans_15 <- int_15$integrate(f, a, b)
  err_15 <- int_15$last_error
  expect_equal(ans_15, area_R, tolerance=eps)

  intervals <- int_15$last_intervals
  tmp <- cbind_list(intervals)
  tmp <- tmp[order(tmp[,1]),]

  ## check that the intervals make sense:
  expect_identical(tmp[1,1], a)
  expect_identical(tmp[nrow(tmp),2], b)
  expect_identical(tmp[-1,1], tmp[-nrow(tmp),2])
  expect_true(all(diff(tmp[,1]) > 0))

  ## Reusing the last intervals is easy:
  expect_identical(int_15$integrate_with_last_intervals(f, a, b), ans_15)
  expect_identical(int_15$last_intervals, intervals)

  ## Now, do the integration with a new set of intervals:
  mid <- (a + b)/2
  intervals_2 <- list(c(a, mid), c(mid, b))
  ans_2 <- int_15$integrate_with_intervals(f, intervals_2)
  err_2 <- int_15$last_error
  expect_gt(err_2, err_15)
  expect_false(identical(err_15, err_2))
  expect_equal(int_15$last_intervals, intervals_2)

  ## Then with the computed intervals;
  ans_redo <- int_15$integrate_with_intervals(f, intervals)
  err_redo <- int_15$last_error
  expect_identical(ans_redo, ans_15)
  expect_equal(err_redo, err_15)
  expect_identical(int_15$last_intervals, intervals)

  ## Rescaling intervals is hard to test well, but here goes:
  b1 <- 49.5
  ans_scal <- int_15$integrate_with_last_intervals(f, 0, b1)

  intervals_scal <- int_15$last_intervals
  tmp_scal <- cbind_list(intervals_scal)
  tmp_scal <- tmp_scal[order(tmp_scal[,1]),]
  expect_identical(dim(tmp_scal), dim(tmp))
  expect_equal(tmp_scal, tmp * b1 / b)
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

    expect_true(int_a$is_adaptive)
    expect_false(int_f$is_adaptive)

    ans_a <- int_a$integrate(f, a, b)
    ans_f <- int_f$integrate(f, a, b)
    ans_q <- int_q$integrate(f, a, b)

    expect_identical(ans_f, ans_q)
    expect_identical(int_f$last_error, int_q$last_error)
    if (length(int_a$last_intervals[[1]]) > 1L) {
      expect_false(identical(ans_a, ans_q))
    } else {
      expect_identical(ans_a, ans_q)
    }
  }
})
