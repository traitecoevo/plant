context("Build_schedule")

strategy_types <- get_list_of_strategy_types()
hyperpar_functions <- get_list_of_hyperpar_functions()

test_that("Corner case", {
  for (i in 1:length(strategy_types)) {
    st <- names(strategy_types)[[i]]
    hyperpar <- hyperpar_functions[[i]]
    p <- scm_base_parameters(st)
    expect_error(build_schedule(p), "no residents")
    p <- expand_parameters(trait_matrix(0.1, "lma"), p, hyperpar)
    expect_error(build_schedule(p), "no residents")
  }
})

## TODO: Not yet done.
test_that("Schedule building", {
  for (x in names(strategy_types)[[1]]) {
    ## This is a really dumb test but it should act as a regression test
    ## at least.
    p <- scm_base_parameters(x)
    p$strategies <- list(strategy_types[[x]]())
    p$seed_rain <- 0.1
    p <- build_schedule(p)
    expect_equal(length(p$cohort_schedule_times_default), 141)
    expect_equal(length(p$cohort_schedule_times[[1]]), 176)
  }
})
