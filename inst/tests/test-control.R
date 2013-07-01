source("helper-tree.R")

context("Control")

ctrl <- new(Control)

obj <- ctrl$parameters

expected <- list(
  cohort_gradient_eps = 1e-6,
  cohort_gradient_richardson = as.numeric(FALSE),
  cohort_gradient_richardson_depth = as.numeric(4L),
  environment_light_max_depth= as.numeric(16),
  environment_light_nbase = as.numeric(17),
  environment_light_tol = 1e-6,
  ode_step_size_max = 1e-1,
  ode_step_size_min = 1e-6,
  ode_tol_abs = 1e-6,
  ode_tol_dydt = 0.0,
  ode_tol_rel = 1e-6,
  ode_tol_y = 1.0,
  plant_assimilation_iterations = as.numeric(1000L),
  plant_assimilation_over_distribution = as.numeric(FALSE),
  plant_assimilation_tol = as.numeric(1e-6),
  plant_seed_iterations = as.numeric(1000L),
  plant_seed_tol = as.numeric(1e-6)
  )

keys <- sort(names(expected))
expect_that(sort(names(obj)),
            is_identical_to(sort(names(expected))))
expect_that(obj[keys], is_identical_to(expected[keys]))

## This just checks that we can pull the ode control parameter out.
## Currently we can't actually work with it.
expect_that(ctrl$ode_control, is_a("Rcpp_OdeControl"))
