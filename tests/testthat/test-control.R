context("Control")

test_that("Defaults", {
  expected <- list(
    node_gradient_eps = 1e-6,
    node_gradient_direction = 1L,
    node_gradient_richardson = FALSE,
    node_gradient_richardson_depth = 4, # size_t, so not int
    
    ode_a_dydt = 0.0,
    ode_a_y = 1.0,
    ode_step_size_initial = 1e-6,
    ode_step_size_max = 1e-1,
    ode_step_size_min = 1e-6,
    ode_tol_abs = 1e-6,
    ode_tol_rel = 1e-6,

    function_integration_rule = 21, # size_t so not int

    offspring_production_iterations = 1000, # size_t
    offspring_production_tol = 1e-8, # 1e-6, Had to change this...
    
    save_RK45_cache = FALSE,

    schedule_nsteps   = 20, # size_t
    schedule_eps      = 1e-3,
    schedule_verbose  = FALSE,
    
    newton_tol_abs = 1e-3,
    GSS_tol_abs = 1e-3,
    vulnerability_curve_ncontrol = 1e2,
    ci_abs_tol = 1e-3,
    ci_niter = 1000
  )

  keys <- sort(names(expected))

  ctrl <- Control()
  expect_is(ctrl, "Control")

  expect_identical(sort(names(ctrl)), keys)
  expect_identical(unclass(ctrl)[keys], expected[keys])
})
