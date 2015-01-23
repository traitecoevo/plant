if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("PlantRunner")

test_that("PlantRunner", {
  p <- Plant(Strategy())
  env <- test_environment(10)

  pr <- PlantRunner(p, env)
  expect_that(pr, is_a("PlantRunner"))
  expect_that(pr$plant, is_a("Plant"))
  expect_that(pr$plant$vars_size, is_identical_to(p$vars_size))
})
