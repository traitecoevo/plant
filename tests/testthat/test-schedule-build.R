context("Schedule_build-FF16")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("refine_schedule on no residents runs trivially", {
  for (x in names(strategy_types)) {
    p <- scm_base_parameters(x)
    expect_silent(scm <- run_scm(p, refine_schedule = TRUE))
    expect_equal(length(scm$parameters$strategies), 0)
  }
})

test_that("Schedule building", {
  for (x in c("FF16")) {
    ## This is a really dumb test but it should act as a regression test
    ## at least. Refinement now happens entirely in C++ via
    ## run_scm(refine_schedule = TRUE) -> SCM::refine_schedule().
    p <- scm_base_parameters(x)
    p$strategies <- list(strategy_types[[x]]())
    p$strategies[[1]]$birth_rate_y <- 0.1

    env <- Environment(x)
    ctrl <- Control()

    scm <- run_scm(p, env, ctrl, refine_schedule = TRUE)
    p2 <- scm$parameters

    expect_equal(length(p2$node_schedule_times_default), 141)
    expect_equal(length(p2$node_schedule_times[[1]]), 186)

    ## Refinement leaves the Parameters self-describing (schedule + ode times).
    expect_equal(length(p2$ode_times), length(scm$ode_times))
    expect_true(length(scm$offspring_production) == 1L)
  }
})
