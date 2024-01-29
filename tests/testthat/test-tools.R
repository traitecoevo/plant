context("Tools")

test_that("fixed_environment", {
  env <- fixed_environment("FF16", 0.5)
  expect_is(env, "FF16_Environment")
  expect_equal(env$shading$spline$xy, cbind(c(0, 75, 150), 0.5))
  expect_equal(env$get_environment_at_height(40), 0.5)
})
