context("Individual utilities")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("Default times", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    ## This is the original function, from tree1:
    cmp_cohort_introduction_times <- function(max_time, multiplier=0.2,
                                              min_step_size=1e-5,
                                              max_step_size=2.0) {
      if (min_step_size <= 0)
        stop("The minimum step size must be greater than zero")
      dt <- time <- times <- 0
      while (time <= max_time) {
        dt <- 2^floor(log2(time * multiplier))
        time <- time + max(min(dt, max_step_size), min_step_size)
        times <- c(times, time)
      }
      # Trucate last time to max_time; it may have overshot.
      last(times) <- max_time
      times
    }

    t1 <- 10.123
    tt <- cohort_schedule_times_default(t1)
    expect_identical(tt[[1]], 0.0)

    expect_lt(last(tt), t1)
    expect_equal(c(tt, t1), cmp_cohort_introduction_times(t1))
  }
})

test_that("Cohort schedule max time", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    t <- cohort_schedule_max_time_default(p)
    d <- Disturbance(p$disturbance_mean_interval)
    expect_equal(t, d$cdf(p$control$schedule_patch_survival))
  }
})

test_that("Default schedule", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]](), strategy_types[[x]]()),
      seed_rain=c(pi/2, pi),
      is_resident=c(TRUE, TRUE))
    cohort_schedule <- cohort_schedule_default(p)
    expect_is(cohort_schedule, "CohortSchedule")
    expect_equal(cohort_schedule$n_species, length(p$strategies))
    t_max <- cohort_schedule_max_time_default(p)
    tt <- cohort_schedule_times_default(t_max)
    expect_identical(cohort_schedule$times(1), tt)
    expect_identical(cohort_schedule$times(2), tt)
  }
})

test_that("strategy_list", {
  for (x in c("FF16", "FF16r")) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    s <- strategy_list(trait_matrix(1, "lma"), p, make_hyperpar(x)())
    expect_equal(length(s), 1)
    expect_is(s, "list")
    expect_is(s[[1]], sprintf("%s_Strategy", x))
  }
})

test_that("individual_list", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()

    obj <- individual_list(trait_matrix(1, "lma"), p, make_hyperpar(x)())
    expect_equal(length(obj), 1)
    expect_is(obj, "list")
    expect_is(obj[[1]], "Individual")
    expect_is(obj[[1]], sprintf("Individual<%s,%s>", x, e))
  }
})
