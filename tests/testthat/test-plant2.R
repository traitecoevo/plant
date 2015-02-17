context("Plant2")

test_that("Reference comparison", {
  s <- Strategy()
  p1 <- Plant(s)
  p2 <- Plant2(s)

  expect_that(p2$strategy, is_identical_to(s))
  expect_that(p2, is_a("Plant2"))

  h0 <- 10
  p1$height <- h0
  p2$height <- h0

  vars1 <- p1$internals
  vars2 <- p2$internals

  expect_that(vars2, is_a("Plant2_internals"))

  expect_that(vars2[["height"]], is_identical_to(vars1[["height"]]))
  ## TODO: @dfalster continue adding tests along here, especially for
  ## the physiological calculations.
})
