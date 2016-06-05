context("Uniroot")

test_that("Agrees with R", {
  quadratic_roots <- function(a, b, c) {
    (-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)
  }
  f <- function(x) x*x - 3*x + 1

  sol <- quadratic_roots(1, -3, 1)
  cmp <- c(uniroot(f, c(0, 1))$root,
           uniroot(f, c(1, 3))$root)

  ## Check the R solutions:
  expect_equal(cmp, sol, tolerance=1e-5)

  ## Then with our thing:
  expect_equal(test_uniroot(f, 0, 1), sol[[1]], tolerance=1e-6)
  expect_equal(test_uniroot(f, 1, 3), sol[[2]], tolerance=1e-6)
})
