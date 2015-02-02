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
  p <- Parameters(strategies=list(Strategy()))
  expect_that(p$seed_rain, equals(1.0))
  expect_that(p$is_resident, is_true())

  expect_that(Parameters(seed_rain=pi),
              throws_error("Incorrect length seed_rain"))
  expect_that(Parameters(is_resident=FALSE),
              throws_error("Incorrect length is_resident"))

  expect_that(Parameters(strategies=list(Strategy(), Strategy()),
                         seed_rain=pi),
              throws_error("Incorrect length seed_rain"))
  expect_that(Parameters(strategies=list(Strategy(), Strategy()),
                         is_resident=TRUE),
              throws_error("Incorrect length is_resident"))

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

test_that("Parameters overwrites Strategy control", {
  ctrl <- ctrl_s <- ctrl_p <- Control()
  ## set these just as markers:
  ctrl_s$schedule_eps <- 1
  ctrl_p$schedule_eps <- 2

  s <- Strategy(control=ctrl_s)
  expect_that(s$control, is_identical_to(ctrl_s))
  expect_that(s$control, not(is_identical_to(ctrl_p)))

  p <- Parameters(control=ctrl_p)
  expect_that(p$control, not(is_identical_to(ctrl_s)))
  expect_that(p$control, is_identical_to(ctrl_p))

  p$strategies <- list(s)
  p$seed_rain <- 1
  p$is_resident <- TRUE
  ## Pass though to force validation:
  tmp <- Patch(p)$parameters
  expect_that(tmp$control, is_identical_to(ctrl_p))
  expect_that(tmp$strategies[[1]]$control, is_identical_to(ctrl_p))

  ## In one shot:
  p2 <- Parameters(control=ctrl_p, strategies=list(s), seed_rain=1,
                   is_resident=TRUE)
  expect_that(p2$control, is_identical_to(ctrl_p))
  expect_that(p2$strategies[[1]]$control, is_identical_to(ctrl_p))
})
