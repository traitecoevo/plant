if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("CohortSchedule")

test_that("CohortScheduleEvent", {
  e <- CohortScheduleEvent(pi, 1)

  expect_that(e$species_index, is_identical_to(1L))
  expect_that(e$species_index_raw, is_identical_to(0.0))

  e$species_index <- 2L
  expect_that(e$species_index, is_identical_to(2L))
  expect_that(e$species_index_raw, is_identical_to(1.0))

  expect_that(e$times,             is_identical_to(pi))
  expect_that(e$time_introduction, is_identical_to(pi))
  expect_that(e$time_end,          is_identical_to(pi))
})
