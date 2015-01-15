if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("Cohort")

test_that("misc", {
  s <- Strategy()

  plant <- Plant(s)
  cohort <- Cohort(s)

  expect_that(cohort, is_a("Cohort"))
  expect_that(cohort$plant, is_a("Plant"))
})
