context("Control")

test_that("Defaults", {
  expected <- list(
    cohort_gradient_eps = 1e-6,
    cohort_gradient_direction = 1L,
    cohort_gradient_richardson = FALSE,
    cohort_gradient_richardson_depth = 4, # size_t, so not int
    environment_light_max_depth= 16, # size_t
    environment_light_nbase = 17, # size_t
    environment_light_tol = 1e-6,
    environment_rescale_usually = FALSE,
    ode_a_dydt = 0.0,
    ode_a_y = 1.0,
    ode_step_size_initial = 1e-6,
    ode_step_size_max = 1e-1,
    ode_step_size_min = 1e-6,
    ode_tol_abs = 1e-6,
    ode_tol_rel = 1e-6,
    plant_assimilation_adaptive = TRUE,
    plant_assimilation_iterations = 1000, # size_t so not int
    plant_assimilation_rule = 21, # size_t so not int
    plant_assimilation_over_distribution = FALSE,
    plant_assimilation_tol = 1e-6,
    plant_seed_iterations = 1000, # size_t
    plant_seed_tol = 1e-8, # 1e-6, Had to change this...

    schedule_nsteps   = 20, # size_t
    schedule_eps      = 1e-3,
    schedule_verbose  = FALSE,
    schedule_patch_survival = 6.25302620663814e-05,

    equilibrium_nsteps   = 20, # size_t
    equilibrium_eps      = 1e-5,
    equilibrium_large_seed_rain_change = 10.0,
    equilibrium_verbose  = TRUE,
    equilibrium_solver_name = "iteration",
    equilibrium_extinct_seed_rain = 1e-3,
    equilibrium_nattempts   = 5, # size_t
    equilibrium_solver_logN = TRUE,
    equilibrium_solver_try_keep = TRUE)

  keys <- sort(names(expected))

  ctrl <- Control()
  expect_is(ctrl, "Control")

  expect_identical(sort(names(ctrl)), keys)
  expect_identical(unclass(ctrl)[keys], expected[keys])
})
