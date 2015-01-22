#include <tree2/control.h>

namespace tree2 {

Control::Control() : integrator(15, 1, 0, 0) {
  plant_assimilation_adaptive = true;

  plant_assimilation_over_distribution = false;
  plant_assimilation_tol = 1e-6;
  plant_assimilation_iterations = 1000;
  plant_assimilation_rule = 21;
  plant_assimilation_reuse_intervals = true;

  // all of this can be deleted?
  plant_assimilation_approximate_use = false;
  plant_assimilation_approximate_tol = 1e-6;
  plant_assimilation_approximate_nbase = 17;
  plant_assimilation_approximate_max_depth = 16;
  plant_assimilation_approximate_akima = false;
  plant_assimilation_approximate_linear = false;
  plant_assimilation_approximate_rescale_usually = false;

  plant_seed_tol = 1e-8;
  plant_seed_iterations = 1000;

  cohort_gradient_eps = 1e-6;
  cohort_gradient_direction = 1;
  cohort_gradient_richardson = false;
  cohort_gradient_richardson_depth = 4;

  environment_light_tol = 1e-6;
  environment_light_nbase = 17;
  environment_light_max_depth = 16;
  environment_light_akima = false; // can be deleted?
  environment_light_linear = false; // can be deleted?
  environment_light_rescale_usually = false;
  environment_light_skip = false; // can be deleted?

  ode_step_size_initial = 1e-6;
  ode_step_size_min = 1e-6;
  ode_step_size_max = 1e-1;
  ode_tol_rel       = 1e-6;
  ode_tol_abs       = 1e-6;
  ode_a_y           = 1.0;
  ode_a_dydt        = 0.0;

  schedule_nsteps   = 20;
  schedule_eps      = 1e-3;
  schedule_progress = false;
  schedule_verbose  = false;
  // This odd number is designed to agree with Daniel's implementation
  // of the model.
  schedule_patch_survival = 6.25302620663814e-05;

  equilibrium_nsteps   = 20;
  equilibrium_eps      = 1e-5;
  equilibrium_large_seed_rain_change = 10;
  equilibrium_progress = false;
  equilibrium_verbose  = true;
  equilibrium_solver   = 1; // TODO: This is a hack
  equilibrium_extinct_seed_rain = 1e-3;
  equilibrium_runsteady_tol = 1e-2;
  equilibrium_inviable_test_eps = 1e-2;
  equilibrium_nattempts = 5;
  equilibrium_solver_logN = true;
  equilibrium_solver_try_keep = true;
}

void Control::initialize() {
  if (!plant_assimilation_adaptive) {
    plant_assimilation_iterations = 1;
  }
  integrator = quadrature::QAG(plant_assimilation_rule,
                               plant_assimilation_iterations,
                               plant_assimilation_tol,
                               plant_assimilation_tol);
}

}
