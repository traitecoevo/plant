context("fitness support")

test_that("positive_1d", {
  f <- function(x) -x^2 + 1
  tol <- 1e-8
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_lt(r[[1]], -1 + tol)
  expect_gt(r[[2]], 1 - tol)

  tol <- 1e-3
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_lt(r[[1]], -1 + tol)
  expect_gt(r[[2]], 1 - tol)

  expect_error(positive_1d(f, -2, 0.1, tol=tol), "no positive values")
})

test_that("positive_2d", {
  skip_if_no_plant_ml_python()
  f <- function(x) {
    if (!is.matrix(x)) {
      x <- rbind(x, deparse.level=0)
    }
    -rowSums(x ^ 2) + 1
  }

  ans <- positive_2d(f, c(0, 0), -2, 2)
})
