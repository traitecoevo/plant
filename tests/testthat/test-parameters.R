if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("Parameters")

test_that("Creation & defaults", {
  p <- Parameters()
  expect_that(p, is_a("Parameters"))

  expect_that(length(p$strategies), equals(0))
  expect_that(length(p$is_resident), equals(0))
  expect_that(length(p$seed_rain), equals(0))

  expected <- list(Pi_0=0.25,
                   c_ext=0.5,
                   n_patches=1,    # NOTE: Different to tree 0.1
                   patch_area=1.0) # NOTE: Different to tree 0.1

  expect_that(p[names(expected)], is_identical_to(expected))
  expect_that(p$disturbance$mean_interval, is_identical_to(30.0))

  expect_that(p$strategy_default, equals(Strategy()))
})
