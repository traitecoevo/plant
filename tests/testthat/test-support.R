context("Support")

test_that("Control presets", {
  ## control() is a lowercase alias for the Control() constructor
  expect_is(control(), "Control")
  expect_equal(control(), Control())

  ## Defaults are the pragmatic, fast-ish settings used for most runs
  expect_equal(Control()$ode_tol_rel, 1e-4)
  expect_equal(Control()$ode_tol_abs, 1e-4)
  expect_equal(Control()$ode_step_size_max, 5)
  expect_equal(Control()$node_gradient_direction, -1)
  expect_equal(Control()$schedule_eps, 2e-2)

  ## control_accurate() tightens the ODE and schedule tolerances
  acc <- control_accurate()
  expect_is(acc, "Control")
  expect_equal(acc$ode_tol_rel, 1e-6)
  expect_equal(acc$ode_tol_abs, 1e-6)
  expect_equal(acc$ode_step_size_max, 1e-1)
  expect_equal(acc$schedule_eps, 1e-3)
})
