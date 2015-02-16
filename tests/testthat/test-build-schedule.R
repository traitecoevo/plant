context("build_schedule")

test_that("Corner case", {
  p <- ebt_base_parameters()
  expect_that(build_schedule(p), throws_error("no residents"))
  p <- expand_parameters(trait_matrix(0.1, "lma"), p)
  expect_that(build_schedule(p), throws_error("no residents"))
})
