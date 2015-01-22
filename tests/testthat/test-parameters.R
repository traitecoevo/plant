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
                   patch_area=1.0, # NOTE: Different to tree 0.1
                   disturbance_mean_interval=30.0)

  expect_that(p[names(expected)], is_identical_to(expected))
  expect_that(p$strategy_default, equals(Strategy()))

  expect_that(p$cohort_schedule_max_time,
              is_identical_to(cohort_schedule_max_time_default(p)))
  expect_that(p$cohort_schedule_times_default,
              is_identical_to(cohort_schedule_times_default(p$cohort_schedule_max_time)))
  expect_that(p$cohort_schedule_times, is_identical_to(list()))
})

test_that("Nontrivial creation", {
  expect_that(Parameters(strategies=list(Strategy())),
              throws_error("Inconsistent lengths"))
  expect_that(Parameters(strategies=list(Strategy()),
                         seed_rain=pi),
              throws_error("Inconsistent lengths"))
  p <- Parameters(strategies=list(Strategy()),
                  seed_rain=pi,
                  is_resident=TRUE)

  expect_that(p$cohort_schedule_max_time,
              is_identical_to(cohort_schedule_max_time_default(p)))
  expect_that(p$cohort_schedule_times_default,
              is_identical_to(cohort_schedule_times_default(p$cohort_schedule_max_time)))
  expect_that(p$cohort_schedule_times,
              is_identical_to(list(p$cohort_schedule_times_default)))

  ## Now, with some of these set:
  t1 <- 10.123
  p <- Parameters(strategies=list(Strategy()),
                  seed_rain=pi,
                  is_resident=TRUE,
                  cohort_schedule_max_time=t1)

  expect_that(p$cohort_schedule_max_time,
              is_identical_to(t1))
  expect_that(p$cohort_schedule_times_default,
              is_identical_to(cohort_schedule_times_default(t1)))
  expect_that(p$cohort_schedule_times,
              is_identical_to(list(p$cohort_schedule_times_default)))

  tt <- p$cohort_schedule_times_default
  p <- Parameters(strategies=list(Strategy()),
                  seed_rain=pi,
                  is_resident=TRUE,
                  cohort_schedule_times_default=tt)

  ## This might be an alarming gap here:
  expect_that(p$cohort_schedule_max_time,
              is_identical_to(cohort_schedule_max_time_default(p)))
  expect_that(p$cohort_schedule_times_default,
              is_identical_to(tt))
  expect_that(p$cohort_schedule_times,
              is_identical_to(list(tt)))
})
