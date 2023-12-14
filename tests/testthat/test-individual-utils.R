context("Individual utilities")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("Default times", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    ## This is the original function, from tree1:
    cmp_node_introduction_times <- function(max_time, multiplier=0.2,
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
    tt <- node_schedule_times_default(t1)
    expect_identical(tt[[1]], 0.0)

    expect_lt(last(tt), t1)
    expect_equal(c(tt, t1), cmp_node_introduction_times(t1))
  }
})

test_that("Node schedule max time", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    expect_lt(max(p$node_schedule_times_default), p$max_patch_lifetime)
  }
})

test_that("Default schedule", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]](), strategy_types[[x]]()))
    node_schedule <- node_schedule_default(p)
    expect_is(node_schedule, "NodeSchedule")
    expect_equal(node_schedule$n_species, length(p$strategies))
    tt <- node_schedule_times_default(p$max_patch_lifetime)
    expect_identical(node_schedule$times(1), tt)
    expect_identical(node_schedule$times(2), tt)
  }
})

test_that("strategy_list", {
  for (x in c("FF16", "FF16r", "FF16drivers")) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    s <- strategy_list(trait_matrix(1, "lma"), p, make_hyperpar(x)(), 1.0)
    expect_equal(length(s), 1)
    expect_is(s, "list")
    expect_is(s[[1]], sprintf("%s_Strategy", x))
  }
})

test_that("individual_list", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    
    obj <- individual_list(trait_matrix(1, "lma"), p, make_hyperpar(x)(), 1.0)
    expect_equal(length(obj), 1)
    expect_is(obj, "list")
    expect_is(obj[[1]], "Individual")
    expect_is(obj[[1]], sprintf("Individual<%s,%s>", x, e))
  }
})
