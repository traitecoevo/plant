source("helper-tree.R")

context("Control")

ctrl <- new(Control)

obj <- ctrl$parameters

expected <- list(
  cohort_gradient_eps = 1e-6,
  cohort_gradient_direction = as.numeric(1L),
  cohort_gradient_richardson = as.numeric(FALSE),
  cohort_gradient_richardson_depth = as.numeric(4L),
  environment_light_max_depth= as.numeric(16),
  environment_light_nbase = as.numeric(17),
  environment_light_tol = 1e-6,
  environment_light_akima = as.numeric(FALSE),
  environment_light_linear = as.numeric(FALSE),
  environment_light_rescale_usually = as.numeric(FALSE),
  environment_light_skip = as.numeric(FALSE),
  ode_a_dydt = 0.0,
  ode_a_y = 1.0,
  ode_step_size_max = 1e-1,
  ode_step_size_min = 1e-6,
  ode_tol_abs = 1e-6,
  ode_tol_rel = 1e-6,
  plant_assimilation_adaptive = as.numeric(TRUE),
  plant_assimilation_iterations = as.numeric(1000L),
  plant_assimilation_rule = as.numeric(21),
  plant_assimilation_over_distribution = as.numeric(FALSE),
  plant_assimilation_reuse_intervals = as.numeric(TRUE),
  plant_assimilation_tol = as.numeric(1e-6),
  plant_assimilation_approximate_use = as.numeric(FALSE),
  plant_assimilation_approximate_max_depth= as.numeric(16),
  plant_assimilation_approximate_nbase = as.numeric(17),
  plant_assimilation_approximate_tol = 1e-6,
  plant_assimilation_approximate_akima = as.numeric(FALSE),
  plant_assimilation_approximate_linear = as.numeric(FALSE),
  plant_assimilation_approximate_rescale_usually = as.numeric(FALSE),
  plant_seed_iterations = as.numeric(1000L),
  plant_seed_tol = as.numeric(1e-6),

  schedule_nsteps   = as.numeric(20L),
  schedule_eps      = 1e-3,
  schedule_progress = as.numeric(FALSE),
  schedule_verbose  = as.numeric(FALSE),
  schedule_default_patch_survival = 6.25302620663814e-05,
  schedule_default_multipler     = 0.2,
  schedule_default_min_step_size = 1e-5,
  schedule_default_max_step_size = 2.0,

  equilibrium_nsteps   = as.numeric(20L),
  equilibrium_eps      = 1e-5,
  equilibrium_large_seed_rain_change = 10.0,
  equilibrium_progress = as.numeric(FALSE),
  equilibrium_verbose  = as.numeric(TRUE),
  equilibrium_solver  = as.numeric(1L),
  equilibrium_extinct_seed_rain = 1e-3,
  equilibrium_runsteady_tol = 1e-2,
  equilibrium_inviable_test_eps = 1e-2,
  equilibrium_nattempts   = as.numeric(5L),
  equilibrium_solver_logN = as.numeric(TRUE),
  equilibrium_solver_try_keep = as.numeric(TRUE)
  )

keys <- sort(names(expected))
test_that("Control keys are as expected", {
  expect_that(sort(names(obj)), is_identical_to(keys))
  expect_that(obj[keys], is_identical_to(expected[keys]))
})

## This just checks that we can pull the ode control parameter out.
test_that("We can pull ode_control out correctly", {
  expect_that(ctrl$ode_control, is_a("Rcpp_OdeControl"))
})

## Check that the ODE control object has the expected parameters
expected.ode <- expected[grep("^ode_", names(expected), value=TRUE)]
names(expected.ode) <- sub("^ode_", "", names(expected.ode))
keys.ode <- sort(names(expected.ode))

obj.ode <- ctrl$ode_control$parameters
test_that("ode_control's keys are as expected", {
  expect_that(sort(names(obj.ode)),
              is_identical_to(sort(names(expected.ode))))
  expect_that(obj.ode[keys.ode], is_identical_to(expected.ode[keys.ode]))
})
