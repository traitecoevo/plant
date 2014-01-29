source("helper-tree.R")

context("Quadrature")

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
  int.15 <- new(QK, 15)
  expect_that(int.15$integrate(g, a, b), equals(area.R))
  expect_that(int.15$last_error,         equals(error.R))
  int.21 <- new(QK, 21)
  expect_that(int.21$integrate(g, a, b), equals(area.R))
  expect_that(int.21$last_error,         equals(error.R))
  int.31 <- new(QK, 31)
  expect_that(int.31$integrate(g, a, b), equals(area.R))
  expect_that(int.31$last_error,         equals(error.R))
  int.41 <- new(QK, 41)
  expect_that(int.41$integrate(g, a, b), equals(area.R))
  expect_that(int.41$last_error,         equals(error.R))
  int.51 <- new(QK, 51)
  expect_that(int.51$integrate(g, a, b), equals(area.R))
  expect_that(int.51$last_error,         equals(error.R))
  int.61 <- new(QK, 61)
  expect_that(int.61$integrate(g, a, b), equals(area.R))
  expect_that(int.61$last_error,         equals(error.R))

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
  expect_that(new(QK, 0),   throws_error())
  expect_that(new(QK, 100), throws_error())
  expect_that(new(QK, NA),  throws_error())
})

test_that("Vectorised interface to integration works", {
  int.15 <- new(QK, 15)
  x <- int.15$integrate_vector_x(a, b)
  expect_that(int.15$integrate_vector(f(x), a, b),
              is_identical_to(int.15$integrate(g, a, b)))

  int.21 <- new(QK, 21)
  x <- int.21$integrate_vector_x(a, b)
  expect_that(int.21$integrate_vector(f(x), a, b),
              is_identical_to(int.21$integrate(g, a, b)))

  int.31 <- new(QK, 31)
  x <- int.31$integrate_vector_x(a, b)
  expect_that(int.31$integrate_vector(f(x), a, b),
              is_identical_to(int.31$integrate(g, a, b)))

  int.41 <- new(QK, 41)
  x <- int.41$integrate_vector_x(a, b)
  expect_that(int.41$integrate_vector(f(x), a, b),
              is_identical_to(int.41$integrate(g, a, b)))

  int.51 <- new(QK, 51)
  x <- int.51$integrate_vector_x(a, b)
  expect_that(int.51$integrate_vector(f(x), a, b),
              is_identical_to(int.51$integrate(g, a, b)))

  int.61 <- new(QK, 61)
  x <- int.61$integrate_vector_x(a, b)
  expect_that(int.61$integrate_vector(f(x), a, b),
              is_identical_to(int.61$integrate(g, a, b)))

  ## Safe from wrong-length arguments.
  expect_that(int.61$integrate_vector(f(x)[-1], a, b),
              throws_error())
  expect_that(int.61$integrate_vector(c(1, f(x)), a, b),
              throws_error())
})
