context("OdeControl")

test_that("Defaults", {
  expected <- list(
    tol_abs=1e-8,
    tol_rel=1e-8,
    a_y=1.0,
    a_dydt=0.0,
    step_size_min=1e-8,
    step_size_max=10.0,
    step_size_initial=1e-6)
  keys <- sort(names(expected))

  ctrl <- OdeControl()
  expect_is(ctrl, "OdeControl")

  expect_identical(sort(names(ctrl)), keys)
  expect_identical(unclass(ctrl)[keys], expected[keys])
})
