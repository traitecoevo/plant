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
                   disturbance_mean_interval=30.0,
                   hyperpar=NULL)

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
  p <- Parameters(strategies=list(Strategy()),
                  seed_rain=pi,
                  is_resident=TRUE,
                  disturbance_mean_interval=2)

  expect_that(p$cohort_schedule_max_time,
              is_less_than(10))
  expect_that(p$cohort_schedule_times_default,
              is_identical_to(cohort_schedule_times_default(p$cohort_schedule_max_time)))
  expect_that(p$cohort_schedule_times,
              is_identical_to(list(p$cohort_schedule_times_default)))
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

test_that("Generate cohort schedule", {
  p <- Parameters(strategies=list(Strategy()),
                  seed_rain=pi/2,
                  is_resident=TRUE)
  sched <- make_cohort_schedule(p)

  expect_that(sched$n_species, equals(1))
  expect_that(sched$max_time, equals(p$cohort_schedule_max_time))
  expect_that(sched$times(1), equals(p$cohort_schedule_times_default))
  expect_that(sched$all_times, equals(p$cohort_schedule_times))
})

test_that("Store hyperparams", {
  p <- Parameters(hyperpar=ff_parameters)
  expect_that(p$hyperpar, is_identical_to(ff_parameters))
  tmp <- Patch(p)$parameters
  expect_that(tmp$hyperpar, is_identical_to(ff_parameters))
})

test_that("ebt_base_parameters", {
  p <- ebt_base_parameters()
  expect_that(p$hyperpar, equals(ff_parameters))
})

test_that("Disturbance interval", {
  p <- ebt_base_parameters()
  expect_that(p$disturbance_mean_interval, equals(30.0))
  expect_that(p$cohort_schedule_max_time,
              is_identical_to(cohort_schedule_max_time_default(p)))
  p$strategies <- list(Strategy())

  p$disturbance_mean_interval <- 10.0
  ## This is going to force us back through the validator
  p2 <- validate(p)
  expect_that(p2$cohort_schedule_max_time,
              is_identical_to(cohort_schedule_max_time_default(p2)))
  expect_that(p2$cohort_schedule_max_time,
              is_less_than(p$cohort_schedule_max_time))
  expect_that(last(p2$cohort_schedule_times_default),
              is_less_than(p2$cohort_schedule_max_time))
  expect_that(last(p2$cohort_schedule_times_default),
              is_less_than(p2$cohort_schedule_max_time))
  expect_that(p2$cohort_schedule_times,
              equals(list(p2$cohort_schedule_times_default)))

  ## We will blow away any data that is stored in p$cohort_schedule*
  expect_that(p$cohort_schedule_max_time, is_more_than(100))
  p$cohort_schedule_max_time <- 1
  p$cohort_schedule_times_default <- 1:10
  p$cohort_schedule_time <- list(1:11)
  expect_that(validate(p), equals(p2))
})
