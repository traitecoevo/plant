context("Schedule_build-FF16")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("Corner case", {
  for (x in names(strategy_types)) {
    p <- scm_base_parameters(x)
    expect_error(build_schedule(p), "no residents")
    p <- expand_parameters(trait_matrix(0.1, "lma"), p)
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
    p$birth_rate <- 0.1
    env <- make_environment(x)
    ctrl <- scm_base_control()

    res <- build_schedule(p, env, ctrl) # state = NULL since no arg provided

    # This has changed, perhaps due split_times missing a cohort
    pars <- res$parameters
    expect_equal(length(pars$cohort_schedule_times_default), 141)
    expect_equal(length(pars$cohort_schedule_times[[1]]), 176)

    # not really relevant, just showing it's here
    state <- res$state
    expect_null(state)
  }
})
