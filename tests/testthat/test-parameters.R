context("Parameters")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("Creation & defaults", {
  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    expect_is(p, sprintf("Parameters<%s,%s>", x, e))

    expect_equal(length(p$strategies), 0)
    expect_equal(length(p$is_resident), 0)
    expect_equal(length(p$seed_rain), 0)

    expected <- list(k_I=0.5,
                     n_patches=1,    # NOTE: Different to tree 0.1
                     patch_area=1.0, # NOTE: Different to tree 0.1
                     disturbance_mean_interval=30.0)

    expect_equal(p[names(expected)], expected)
    expect_equal(p$strategy_default, s)

    expect_identical(p$cohort_schedule_max_time, cohort_schedule_max_time_default(p))
    expect_identical(p$cohort_schedule_times_default, cohort_schedule_times_default(p$cohort_schedule_max_time))
    expect_identical(p$cohort_schedule_times, list())
  }
})

test_that("Nontrivial creation", {
  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(s))
    expect_equal(p$seed_rain, 1.0)
    expect_true(p$is_resident)

    expect_error(Parameters(x, e)(seed_rain=pi), "Incorrect length seed_rain")
    expect_error(Parameters(x, e)(is_resident=FALSE), "Incorrect length is_resident")

    expect_error(Parameters(x, e)(strategies=list(strategy_types[[x]](),
                                   strategy_types[[x]]()),
                                 seed_rain=pi),
                 "Incorrect length seed_rain")
    expect_error(Parameters(x, e)(strategies=list(strategy_types[[x]](),
                                   strategy_types[[x]]()),
                                 is_resident=TRUE),
                 "Incorrect length is_resident")

    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi,
                          is_resident=TRUE)

    expect_identical(p$cohort_schedule_max_time, cohort_schedule_max_time_default(p))
    expect_identical(p$cohort_schedule_times_default, cohort_schedule_times_default(p$cohort_schedule_max_time))
    expect_identical(p$cohort_schedule_times, list(p$cohort_schedule_times_default))

    ## Now, with some of these set:
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi,
                          is_resident=TRUE,
                          disturbance_mean_interval=2)

    expect_lt(p$cohort_schedule_max_time, 10)
    expect_identical(p$cohort_schedule_times_default, cohort_schedule_times_default(p$cohort_schedule_max_time))
    expect_identical(p$cohort_schedule_times, list(p$cohort_schedule_times_default))
  }
})

test_that("Parameters overwrites Strategy control", {
  for (x in names(strategy_types)) {
    ctrl <- ctrl_s <- ctrl_p <- Control()
    ## set these just as markers:
    ctrl_s$schedule_eps <- 1
    ctrl_p$schedule_eps <- 2

    s <- strategy_types[[x]](control=ctrl_s)
    e <- environment_types[[x]]
    expect_identical(s$control, ctrl_s)
    expect_false(identical(s$control, ctrl_p))

    p <- Parameters(x, e)(control=ctrl_p)
    expect_false(identical(p$control, ctrl_s))
    expect_identical(p$control, ctrl_p)

    p$strategies <- list(s)
    p$seed_rain <- 1
    p$is_resident <- TRUE
    ## Pass though to force validation:
    tmp <- Patch(x, e)(p)$parameters
    expect_identical(tmp$control, ctrl_p)
    expect_identical(tmp$strategies[[1]]$control, ctrl_p)

    ## In one shot:
    p2 <- Parameters(x, e)(control=ctrl_p, strategies=list(s),
                           seed_rain=1, is_resident=TRUE)
    expect_identical(p2$control, ctrl_p)
    expect_identical(p2$strategies[[1]]$control, ctrl_p)
  }
})

test_that("Generate cohort schedule", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi/2, is_resident=TRUE)
    sched <- make_cohort_schedule(p)

    expect_equal(sched$n_species, 1)
    expect_equal(sched$max_time, p$cohort_schedule_max_time)
    expect_equal(sched$times(1), p$cohort_schedule_times_default)
    expect_equal(sched$all_times, p$cohort_schedule_times)
  }
})

test_that("Validate", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    expect_equal(validate(p), p)
    p$is_resident <- TRUE
    expect_error(validate(p), "Incorrect length is_resident")
  }
})


test_that("scm_base_parameters", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- scm_base_parameters(x)
    expect_is(p, sprintf("Parameters<%s,%s>", x, e))
  }
})

test_that("Disturbance interval", {
  for (x in names(strategy_types)[[1]]) {
    p <- scm_base_parameters(x)
    expect_equal(p$disturbance_mean_interval, 30.0)
    expect_identical(p$cohort_schedule_max_time, cohort_schedule_max_time_default(p))
    p$strategies <- list(strategy_types[[x]]())

    p$disturbance_mean_interval <- 10.0
    ## This is going to force us back through the validator
    p2 <- validate(p)
    expect_identical(p2$cohort_schedule_max_time, cohort_schedule_max_time_default(p2))
    expect_lt(p2$cohort_schedule_max_time, p$cohort_schedule_max_time)
    expect_lt(last(p2$cohort_schedule_times_default), p2$cohort_schedule_max_time)
    expect_lt(last(p2$cohort_schedule_times_default), p2$cohort_schedule_max_time)
    expect_equal(p2$cohort_schedule_times, list(p2$cohort_schedule_times_default))

    ## We will blow away any data that is stored in p$cohort_schedule*
    expect_gt(p$cohort_schedule_max_time, 100)
    p$cohort_schedule_max_time <- 1
    p$cohort_schedule_times_default <- 1:10
    p$cohort_schedule_time <- list(1:11)
    expect_equal(validate(p), p2)
  }
})
