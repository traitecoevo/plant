context("build_schedule")

test_that("Corner case", {
  p <- ebt_base_parameters()
  expect_that(build_schedule(p), throws_error("no residents"))
  p <- expand_parameters(trait_matrix(0.1, "lma"), p)
  expect_that(build_schedule(p), throws_error("no residents"))
})

test_that("Schedule building", {
  ## This is a really dumb test but it should act as a regression test
  ## at least.
  p <- ebt_base_parameters()
  p$strategies <- list(Strategy())
  p$seed_rain <- 0.1
  p <- build_schedule(p)
  expect_that(length(p$cohort_schedule_times_default), equals(141))
  expect_that(length(p$cohort_schedule_times[[1]]),    equals(176))
})
