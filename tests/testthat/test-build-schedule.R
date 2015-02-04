if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

test_that("Corner case", {
  p <- ebt_base_parameters()
  expect_that(build_schedule(p), throws_error("no residents"))
  p <- expand_parameters(trait_matrix(0.1, "lma"), p)
  expect_that(build_schedule(p), throws_error("no residents"))
})
