context("Trapezium integration")

test_that("Trapezium rule works", {
  n <- 20
  set.seed(1)
  xx <- sort(runif(n))
  yy <- runif(n)
  expect_equal(trapezium(xx, yy), sum((xx[-1] - xx[-n]) * (yy[-1] + yy[-n])) / 2)
  expect_identical(trapezium(c(1, 1), c(1, 2)), 0.0)
})

test_that("Trapezium local error estimate is correct", {
  ## Here is the naive implementation of integration error...
  local_error_integration_R <- function(x, y) {
    n <- length(x)
    i0 <- -c(n-1, n)
    i1 <- -c(1,   n)
    i2 <- -c(1,   2)

    y1.fine <-
      0.5*((y[i0] + y[i1])*(x[i1] - x[i0]) +
           (y[i1] + y[i2])*(x[i2] - x[i1]))
    y1.coarse <-
      0.5*((y[i0] + y[i2])*(x[i2] - x[i0]))
    
    c(NA, abs(y1.coarse - y1.fine), NA)
  }
  
  ## Some data:
  set.seed(1)
  xx <- sort(runif(100, 0, 6*pi))
  yy <- 6*(sin(xx) + 1)
  tot <- abs(trapezium(xx, yy))
  
  expect_equal(local_error_integration(xx, yy, 1), local_error_integration_R(xx, yy))
  expect_equal(local_error_integration(xx, yy, tot), local_error_integration_R(xx, yy)/tot)
  
  expect_equal(local_error_integration(numeric(0), numeric(0), tot), numeric(0))
  expect_equal(local_error_integration(xx[1], yy[1], tot), NA_real_)
  expect_equal(local_error_integration(xx[1:2], yy[1:2], tot), c(NA_real_, NA_real_))
  expect_error(local_error_integration(xx[1:2], yy[1:3], tot))
})
