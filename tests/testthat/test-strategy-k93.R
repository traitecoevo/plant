# Built from  tests/testthat/test-strategy-ff16.R on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16
context("Strategy-K93")

test_that("Defaults", {
  expected <- list(
    height_0 = 2.0,
    size_c = 30.0,
    b_0 = 0.059,
    b_1 = 0.012,
    b_2 = 0.00041,
    c_0 = 0.008,
    c_1 = 0.00044,
    d_0 = 0.00073,
    d_1 = 0.044,
    control = Control())

  keys <- sort(names(expected))

  s <- K93_Strategy()
  expect_is(s, "K93_Strategy")

  expect_identical(sort(names(s)), keys)
  expect_identical(unclass(s)[keys], expected[keys])
})


test_that("Reference comparison", {
  s <- K93_Strategy()
  p <- K93_Plant(s)

  expect_identical(p$strategy, s)

  ## Set the size to something (here 10)
  s0 <- 10
  p$set_state("size", s0)


  expect_identical(p$state("size"), s0)

  ## Check: Is this redundant now
  ## We now use 
  vars <- p$internals
  expect_identical(p$state("size"), vars$states[which(p$ode_names == "size")])
})


test_that("Critical Names", {
  s <- K93_Strategy()
  my_names <- K93_Plant(s)$ode_names
  expect_identical(my_names[1:3], c("size", "mortality", "fecundity"))
})

