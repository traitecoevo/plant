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
    
    env <- make_environment(x)
    ctrl <- scm_base_control()

    p <- build_schedule(p, env, ctrl)
    expect_equal(length(p$node_schedule_times_default), 141)
    expect_equal(length(p$node_schedule_times[[1]]), 176)
  }
})
