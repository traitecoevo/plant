context("Tools")
strategy_types <- get_list_of_strategy_types()

test_that("fixed_environment", {
  env <- fixed_environment(0.5)
  expect_is(env, "LightEnvironment")
    expect_equal(env$environment_interpolator$xy, cbind(c(0, 75, 150), 0.5))
  expect_equal(env$canopy_openness(40), 0.5)
})
