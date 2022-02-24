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

    assimilator_adaptive_integration = TRUE,
    assimilator_integration_iterations = 1000, # size_t so not int
    assimilator_integration_rule = 21, # size_t so not int
    assimilator_integration_tol = 1e-6,

    offspring_production_iterations = 1000, # size_t
    offspring_production_tol = 1e-8, # 1e-6, Had to change this...

    schedule_nsteps   = 20, # size_t
    schedule_eps      = 1e-3,
    schedule_verbose  = FALSE,

    equilibrium_nsteps   = 20, # size_t
    equilibrium_eps      = 1e-5,
    equilibrium_large_birth_rate_change = 10.0,
    equilibrium_verbose  = TRUE,
    equilibrium_solver_name = "iteration",
    equilibrium_extinct_birth_rate = 1e-3,
    equilibrium_nattempts   = 5, # size_t
    equilibrium_solver_logN = TRUE,
    equilibrium_solver_try_keep = TRUE)

  keys <- sort(names(expected))

  ctrl <- Control()
  expect_is(ctrl, "Control")

  expect_identical(sort(names(ctrl)), keys)
  expect_identical(unclass(ctrl)[keys], expected[keys])
})
