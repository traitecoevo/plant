context("Schedule_build-FF16")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("Corner case", {
  for (x in names(strategy_types)) {
    p <- scm_base_parameters(x)
    expect_error(build_schedule(p), "no residents")
  }
})

## TODO: Not yet done.
test_that("Schedule building", {
  for (x in c("FF16")) {
    ## This is a really dumb test but it should act as a regression test
    ## at least.
    p <- scm_base_parameters(x)
    p$strategies <- list(strategy_types[[x]]())
    p$strategies[[1]]$birth_rate_y <- 0.1
    
    env <- Environment(x)
    ctrl <- scm_base_control()

    p <- build_schedule(p, env, ctrl)
    expect_equal(length(p$node_schedule_times_default), 141)
    expect_equal(length(p$node_schedule_times[[1]]), 186)
  }
})

test_that("C++ refine_schedule matches R build_schedule", {
  for (x in c("FF16")) {
    e <- environment_types[[x]]

    p <- scm_base_parameters(x)
    p$strategies <- list(strategy_types[[x]]())
    p$strategies[[1]]$birth_rate_y <- 0.1

    env <- Environment(x)
    ctrl <- scm_base_control()

    p_old <- build_schedule(p, env, ctrl)

    scm <- SCM(x, e)(p, env, ctrl)
    scm$refine_schedule()
    p_new <- scm$parameters

    ## Same refined schedule, ode times, and offspring production, computed
    ## entirely in C++.
    expect_equal(p_new$node_schedule_times, p_old$node_schedule_times)
    expect_equal(p_new$ode_times, p_old$ode_times)
    expect_equal(scm$offspring_production,
                 attr(p_old, "offspring_production"))
    expect_equal(length(p_new$node_schedule_times[[1]]), 186)
  }
})
