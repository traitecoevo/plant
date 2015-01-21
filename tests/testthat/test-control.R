if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("Control")

test_that("Defaults", {
  expected <- list(
    cohort_gradient_eps = 1e-6,
    cohort_gradient_direction = 1L,
    cohort_gradient_richardson = FALSE,
    cohort_gradient_richardson_depth = 4, # size_t, so not int
    environment_light_max_depth= 16L,
    environment_light_nbase = 17L,
    environment_light_tol = 1e-6,
    environment_light_akima = FALSE,
    environment_light_linear = FALSE,
    environment_light_rescale_usually = FALSE,
    environment_light_skip = FALSE,
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
    plant_assimilation_reuse_intervals = TRUE,
    plant_assimilation_tol = 1e-6,
    plant_assimilation_approximate_use = FALSE,
    plant_assimilation_approximate_max_depth= 16L,
    plant_assimilation_approximate_nbase = 17L,
    plant_assimilation_approximate_tol = 1e-6,
    plant_assimilation_approximate_akima = FALSE,
    plant_assimilation_approximate_linear = FALSE,
    plant_assimilation_approximate_rescale_usually = FALSE,
    plant_seed_iterations = 1000L,
    plant_seed_tol = 1e-8, # 1e-6, Had to change this...

    schedule_nsteps   = 20L,
    schedule_eps      = 1e-3,
    schedule_progress = FALSE,
    schedule_verbose  = FALSE,
    schedule_default_patch_survival = 6.25302620663814e-05,
    schedule_default_multipler     = 0.2,
    schedule_default_min_step_size = 1e-5,
    schedule_default_max_step_size = 2.0,

    equilibrium_nsteps   = 20L,
    equilibrium_eps      = 1e-5,
    equilibrium_large_seed_rain_change = 10.0,
    equilibrium_progress = FALSE,
    equilibrium_verbose  = TRUE,
    equilibrium_solver  = 1L,
    equilibrium_extinct_seed_rain = 1e-3,
    equilibrium_runsteady_tol = 1e-2,
    equilibrium_inviable_test_eps = 1e-2,
    equilibrium_nattempts   = 5L,
    equilibrium_solver_logN = TRUE,
    equilibrium_solver_try_keep = TRUE)

  keys <- sort(names(expected))

  ctrl <- Control()
  expect_that(ctrl, is_a("Control"))

  expect_that(sort(names(ctrl)), is_identical_to(keys))
  expect_that(unclass(ctrl)[keys], is_identical_to(expected[keys]))
})
