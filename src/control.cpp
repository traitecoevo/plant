#include <plant/control.h>

namespace plant {

Control::Control() {
  plant_assimilation_adaptive = true;

  plant_assimilation_over_distribution = false;
  plant_assimilation_tol = 1e-6;
  plant_assimilation_iterations = 1000;
  plant_assimilation_rule = 21;

  plant_seed_tol = 1e-8;
  plant_seed_iterations = 1000;

  soil_infiltration_rate = 0.0;
  soil_number_of_depths = 0;

  cohort_gradient_eps = 1e-6;
  cohort_gradient_direction = 1;
  cohort_gradient_richardson = false;
  cohort_gradient_richardson_depth = 4;

  environment_light_tol = 1e-6;
  environment_light_nbase = 17;
  environment_light_max_depth = 16;
  environment_rescale_usually = false;

  ode_step_size_initial = 1e-6;
  ode_step_size_min = 1e-6;
  ode_step_size_max = 1e-1;
  ode_tol_rel       = 1e-6;
  ode_tol_abs       = 1e-6;
  ode_a_y           = 1.0;
  ode_a_dydt        = 0.0;

  schedule_nsteps   = 20;
  schedule_eps      = 1e-3;
  schedule_verbose  = false;

  equilibrium_nsteps   = 20;
  equilibrium_eps      = 1e-5;
  equilibrium_large_birth_rate_change = 10;
  equilibrium_verbose  = true;
  equilibrium_solver_name = "iteration";
  equilibrium_extinct_birth_rate = 1e-3;
  equilibrium_nattempts = 5;
  equilibrium_solver_logN = true;
  equilibrium_solver_try_keep = true;
}

}
