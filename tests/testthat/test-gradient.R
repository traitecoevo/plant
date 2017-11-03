context("Gradient")

test_that("Gradients agree", {
  ## Test the simple finite differencing gradient function.
  gradient_fd_forward <- function(f, x, dx) {
    (f(x + dx) - f(x)) / dx
  }
  gradient_fd_centre <- function(f, x, dx) {
    (f(x + dx/2) - f(x - dx/2)) / dx
  }
  gradient_fd_backward <- function(f, x, dx) {
    (f(x - dx) - f(x)) / (-dx)
  }

  f <- function(x) x*x - 3*x + 1

  dx <- 0.001
  x <- 1

  ## Computing f(x)
  expect_identical(test_gradient_fd1(f, x, dx, 1), gradient_fd_forward(f, x, dx))
  expect_identical(test_gradient_fd1(f, x, dx, 0), gradient_fd_centre(f, x, dx))
  expect_identical(test_gradient_fd1(f, x, dx, -1), gradient_fd_backward(f, x, dx))

  ## Providing f(x)
  expect_identical(test_gradient_fd1(f, x, dx, 1, f(x)), gradient_fd_forward(f, x, dx))
  expect_identical(test_gradient_fd1(f, x, dx, 0, f(x)), gradient_fd_centre(f, x, dx))
  expect_identical(test_gradient_fd1(f, x, dx, -1, f(x)), gradient_fd_backward(f, x, dx))


  d <- 1e-6
  r <- 4L
  method_args <- list(d=d, eps=d)
  expect_equal(test_gradient_richardson(f, x, d, r), numDeriv::grad(f, x, method.args=method_args))
})
