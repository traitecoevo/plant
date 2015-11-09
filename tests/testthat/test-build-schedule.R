context("Build_schedule")

strategy_types <- get_list_of_strategy_types()

test_that("Corner case", {
  for (x in names(strategy_types)) {
    p <- scm_base_parameters(x)
    expect_that(build_schedule(p), throws_error("no residents"))
    p <- expand_parameters(trait_matrix(0.1, "lma"), p)
    expect_that(build_schedule(p), throws_error("no residents"))
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
    expect_that(length(p$cohort_schedule_times_default), equals(141))
    expect_that(length(p$cohort_schedule_times[[1]]),    equals(176))
  }
})
