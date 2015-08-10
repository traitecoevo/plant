context("Plant utilities")

## TODO: Remove ["FFW16"] to test this with all types. But first
## requires issue #162 to be resolved
strategy_types <- get_list_of_strategy_types()

test_that("Default times", {
  for (x in names(strategy_types)) {
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
    expect_that(tt[[1]], is_identical_to(0.0))

    expect_that(last(tt), is_less_than(t1))
    expect_that(c(tt, t1),
                equals(cmp_cohort_introduction_times(t1)))
  }
})

test_that("Cohort schedule max time", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)()
    t <- cohort_schedule_max_time_default(p)
    d <- Disturbance(p$disturbance_mean_interval)
    expect_that(t, equals(d$cdf(p$control$schedule_patch_survival)))
  }
})

test_that("Default schedule", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)(strategies=list(strategy_types[[x]](), strategy_types[[x]]()),
      seed_rain=c(pi/2, pi),
      is_resident=c(TRUE, TRUE))
    cohort_schedule <- cohort_schedule_default(p)
    expect_that(cohort_schedule, is_a("CohortSchedule"))
    expect_that(cohort_schedule$n_species, equals(length(p$strategies)))
    t_max <- cohort_schedule_max_time_default(p)
    tt <- cohort_schedule_times_default(t_max)
    expect_that(cohort_schedule$times(1), is_identical_to(tt))
    expect_that(cohort_schedule$times(2), is_identical_to(tt))
  }
})

test_that("strategy_list", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)()
    s <- strategy_list(trait_matrix(1, "lma"), p)
    expect_that(length(s), equals(1))
    expect_that(s, is_a("list"))
    expect_that(s[[1]], is_a(sprintf("%s_Strategy", x)))
  }
})

test_that("plant_list", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)()

    obj <- plant_list(trait_matrix(1, "lma"), p)
    expect_that(length(obj), equals(1))
    expect_that(obj, is_a("list"))
    expect_that(obj[[1]], is_a("Plant"))
    expect_that(obj[[1]], is_a(sprintf("Plant<%s>", x)))
  }
})
