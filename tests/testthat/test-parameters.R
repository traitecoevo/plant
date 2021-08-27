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
    expect_equal(length(p$birth_rate), 0)
    
    expected <- list(n_patches=1,
                     patch_area=1.0,
                     patch_type = "meta-population",
                     max_patch_lifetime=105.32)

    expect_equal(p[names(expected)], expected)
    expect_equal(p$strategy_default, s)

    expect_identical(p$cohort_schedule_times_default, cohort_schedule_times_default(p$max_patch_lifetime))
    expect_identical(p$cohort_schedule_times, list())
  }
})

test_that("Nontrivial creation", {
  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(s))
    expect_equal(p$birth_rate, 1.0)
    expect_true(p$is_resident)
    
    expect_error(Parameters(x, e)(birth_rate=pi), "Incorrect length birth_rate")
    expect_error(Parameters(x, e)(is_resident=FALSE), "Incorrect length is_resident")

    expect_error(Parameters(x, e)(strategies=list(strategy_types[[x]](),
                                   strategy_types[[x]]()),
                                 birth_rate=pi),
                 "Incorrect length birth_rate")
    expect_error(Parameters(x, e)(strategies=list(strategy_types[[x]](),
                                   strategy_types[[x]]()),
                                 is_resident=TRUE),
                 "Incorrect length is_resident")

    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          birth_rate=pi,
                          is_resident=TRUE)

    expect_identical(p$cohort_schedule_times_default, cohort_schedule_times_default(p$max_patch_lifetime))
    expect_identical(p$cohort_schedule_times, list(p$cohort_schedule_times_default))

    ## Now, with some of these set:
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          birth_rate=pi,
                          is_resident=TRUE,
                          max_patch_lifetime=7.021333)

    message(x, " ", p$max_patch_lifetime)
    expect_lt(p$max_patch_lifetime, 10)
    expect_identical(p$cohort_schedule_times_default, cohort_schedule_times_default(p$max_patch_lifetime))
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
    p$birth_rate <- 1
    p$is_resident <- TRUE
    ## Pass though to force validation:
    tmp <- Patch(x, e)(p)$parameters
    expect_identical(tmp$control, ctrl_p)
    expect_identical(tmp$strategies[[1]]$control, ctrl_p)

    ## In one shot:
    p2 <- Parameters(x, e)(control=ctrl_p, strategies=list(s),
                           birth_rate=1, is_resident=TRUE)
    expect_identical(p2$control, ctrl_p)
    expect_identical(p2$strategies[[1]]$control, ctrl_p)
  }
})

test_that("Generate cohort schedule", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          birth_rate=pi/2, is_resident=TRUE)
    sched <- make_cohort_schedule(p)

    expect_equal(sched$n_species, 1)
    expect_equal(sched$max_time, p$max_patch_lifetime)
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


test_that("scm_default_control", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- scm_default_control(x)
    expect_is(p, sprintf("Parameters<%s,%s>", x, e))
  }
})

test_that("Patch runtime", {
  for (x in names(strategy_types)[[1]]) {
    p <- scm_default_control(x)
    expect_equal(p$max_patch_lifetime, 105.32)
    p$strategies <- list(strategy_types[[x]]())

    ## This is going to force us back through the validator
    p$max_patch_lifetime <- 35.10667
    p2 <- validate(p)
    expect_lt(last(p2$cohort_schedule_times_default), p2$max_patch_lifetime)
    expect_equal(p2$cohort_schedule_times, list(p2$cohort_schedule_times_default))

    ## We will blow away any data that is stored in p$cohort_schedule*
    p$cohort_schedule_times_default <- 1:10
    p$cohort_schedule_time <- list(1:11)
    expect_equal(validate(p), p2)
  }
})
