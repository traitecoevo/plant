context("fitness support")

test_that("positive_1d", {
  f <- function(x) -x^2 + 1
  tol <- 1e-8
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_that(r[[1]], is_less_than(-1 + tol))
  expect_that(r[[2]], is_more_than(1 - tol))

  tol <- 1e-3
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_that(r[[1]], is_less_than(-1 + tol))
  expect_that(r[[2]], is_more_than(1 - tol))

  expect_that(positive_1d(f, -2, 0.1, tol=tol),
              throws_error("no positive values"))
})

test_that("positive_2d", {
  skip_if_not_installed("plant.ml")
  f <- function(x) {
    if (!is.matrix(x)) {
      x <- rbind(x, deparse.level=0)
    }
    -rowSums(x ^ 2) + 1
  }

  ans <- positive_2d(f, c(0, 0), -2, 2)
})


