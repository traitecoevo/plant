context("Parameters")

strategy_types <- get_list_of_strategy_types()

test_that("hyperpar creation", {
  for (x in names(strategy_types)) {
    h <- hyperpar(x)
    expect_is(h, "function")
    expect_equal(names(formals(h)), c("m", "s", "filter"))
    expect_equal(make_hyperpar(x)(), h)
  }
})

test_that("Creation & defaults", {
  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    p <- Parameters(x)()
    expect_is(p, sprintf("Parameters<%s>", x))

    expect_equal(length(p$strategies), 0)
    expect_equal(length(p$is_resident), 0)
    expect_equal(length(p$seed_rain), 0)

    expected <- list(k_I=0.5,
                     n_patches=1,    # NOTE: Different to tree 0.1
                     patch_area=1.0, # NOTE: Different to tree 0.1
                     disturbance_mean_interval=30.0,
                     hyperpar=hyperpar(x))

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
    p <- Parameters(x)(strategies=list(s))
    expect_equal(p$seed_rain, 1.0)
    expect_true(p$is_resident)

    expect_error(Parameters(x)(seed_rain=pi), "Incorrect length seed_rain")
    expect_error(Parameters(x)(is_resident=FALSE), "Incorrect length is_resident")

    expect_error(Parameters(x)(strategies=list(strategy_types[[x]](),
                                   strategy_types[[x]]()),
                                 seed_rain=pi),
                 "Incorrect length seed_rain")
    expect_error(Parameters(x)(strategies=list(strategy_types[[x]](),
                                   strategy_types[[x]]()),
                                 is_resident=TRUE),
                 "Incorrect length is_resident")

    p <- Parameters(x)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi,
                          is_resident=TRUE)

    expect_identical(p$cohort_schedule_max_time, cohort_schedule_max_time_default(p))
    expect_identical(p$cohort_schedule_times_default, cohort_schedule_times_default(p$cohort_schedule_max_time))
    expect_identical(p$cohort_schedule_times, list(p$cohort_schedule_times_default))

    ## Now, with some of these set:
    p <- Parameters(x)(strategies=list(strategy_types[[x]]()),
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
    expect_identical(s$control, ctrl_s)
    expect_not_identical(s$control, ctrl_p)

    p <- Parameters(x)(control=ctrl_p)
    expect_not_identical(p$control, ctrl_s)
    expect_identical(p$control, ctrl_p)

    p$strategies <- list(s)
    p$seed_rain <- 1
    p$is_resident <- TRUE
    ## Pass though to force validation:
    tmp <- Patch(x)(p)$parameters
    expect_identical(tmp$control, ctrl_p)
    expect_identical(tmp$strategies[[1]]$control, ctrl_p)

    ## In one shot:
    p2 <- Parameters(x)(control=ctrl_p, strategies=list(s),
                           seed_rain=1, is_resident=TRUE)
    expect_identical(p2$control, ctrl_p)
    expect_identical(p2$strategies[[1]]$control, ctrl_p)
  }
})

test_that("Generate cohort schedule", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)(strategies=list(strategy_types[[x]]()),
                          seed_rain=pi/2, is_resident=TRUE)
    sched <- make_cohort_schedule(p)

    expect_equal(sched$n_species, 1)
    expect_equal(sched$max_time, p$cohort_schedule_max_time)
    expect_equal(sched$times(1), p$cohort_schedule_times_default)
    expect_equal(sched$all_times, p$cohort_schedule_times)
  }
})

test_that("Store hyperparams", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)(hyperpar=hyperpar(x))
    expect_identical(p$hyperpar, hyperpar(x))
    tmp <- Patch(x)(p)$parameters
    expect_identical(tmp$hyperpar, hyperpar(x))
  }
})

test_that("Validate", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)()
    expect_equal(validate(p), p)
    p$is_resident <- TRUE
    expect_error(validate(p), "Incorrect length is_resident")
  }
})


test_that("scm_base_parameters", {
  for (x in names(strategy_types)) {
    p <- scm_base_parameters(x)
    expect_equal(p$hyperpar, hyperpar(x))
    expect_is(p, sprintf("Parameters<%s>", x))
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

test_that("narea calculation", {
  x <- c(1.38, 3.07, 2.94)
  p0 <- FF16_Parameters()
  m <- trait_matrix(x, "hmat")
  expect_not_warning(sl <- strategy_list(m, p0))

  cmp <- lapply(x, function(xi) strategy(trait_matrix(xi, "hmat"), p0))
  expect_equal(sl, cmp)
})
