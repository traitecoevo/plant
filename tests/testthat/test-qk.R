context("QK")

## These will be used frequently:
f <- sin
a <- 0
b <- 1

test_that("R's integration behaves as expected", {
  ans_R <- integrate(f, a, b)
  expect_equal(ans_R$subdivisions, 1)
})

test_that("Integration agrees with R on simple problem", {
  ans_R <- integrate(f, a, b)
  area_R <- ans_R$value
  error_R <- ans_R$abs.error

  int_15 <- QK(15)
  expect_equal(int_15$integrate(f, a, b), area_R)
  expect_equal(int_15$last_error, error_R)
  int_21 <- QK(21)
  expect_equal(int_21$integrate(f, a, b), area_R)
  expect_equal(int_21$last_error, error_R)
  int_31 <- QK(31)
  expect_equal(int_31$integrate(f, a, b), area_R)
  expect_equal(int_31$last_error, error_R)
  int_41 <- QK(41)
  expect_equal(int_41$integrate(f, a, b), area_R)
  expect_equal(int_41$last_error, error_R)
  int_51 <- QK(51)
  expect_equal(int_51$integrate(f, a, b), area_R)
  expect_equal(int_51$last_error, error_R)
  int_61 <- QK(61)
  expect_equal(int_61$integrate(f, a, b), area_R)
  expect_equal(int_61$last_error, error_R)

  ## These should all be different, but the first might not actually
  ## be different.
  ## expect_false(identical(int_15$last_area, int_21$last_area))
  expect_false(identical(int_21$last_area, int_31$last_area))
  expect_false(identical(int_31$last_area, int_41$last_area))
  expect_false(identical(int_41$last_area, int_51$last_area))
  expect_false(identical(int_51$last_area, int_61$last_area))
})

test_that("Cannot make non-existent rules", {
  expect_error(QK(0), "Unknown rule 0")
  expect_error(QK(100), "Unknown rule 100")
  expect_error(QK(NA), "Unknown rule")
})

test_that("Vectorised interface to integration works", {
  int_15 <- QK(15)
  x <- int_15$integrate_vector_x(a, b)
  expect_equal(int_15$integrate_vector(f(x), a, b), int_15$integrate(f, a, b))

  int_21 <- QK(21)
  x <- int_21$integrate_vector_x(a, b)
  expect_equal(int_21$integrate_vector(f(x), a, b), int_21$integrate(f, a, b))

  int_31 <- QK(31)
  x <- int_31$integrate_vector_x(a, b)
  expect_equal(int_31$integrate_vector(f(x), a, b), int_31$integrate(f, a, b))

  int_41 <- QK(41)
  x <- int_41$integrate_vector_x(a, b)
  expect_equal(int_41$integrate_vector(f(x), a, b), int_41$integrate(f, a, b))

  int_51 <- QK(51)
  x <- int_51$integrate_vector_x(a, b)
  expect_equal(int_51$integrate_vector(f(x), a, b), int_51$integrate(f, a, b))

  int_61 <- QK(61)
  x <- int_61$integrate_vector_x(a, b)
  expect_equal(int_61$integrate_vector(f(x), a, b), int_61$integrate(f, a, b))

  ## Safe from wrong-length arguments.
  expect_error(int_61$integrate_vector(f(x)[-1], a, b), "Incorrect length input")
  expect_error(int_61$integrate_vector(c(1, f(x)), a, b), "Incorrect length input")
})
