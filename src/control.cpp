#include <plant/control.h>

namespace plant {

Control::Control() : integrator(15, 1, 0, 0) {

  // patch settings
  patch_area = 1.0;
  patch_max_lifetime = 105.32; // designed to agree with Daniel's implementation

  // integration rules for net_reproduction
  plant_seed_tol = 1e-8;
  plant_seed_iterations = 1000;

  // integration rules for cohorts
  cohort_gradient_eps = 1e-6;
  cohort_gradient_direction = 1;
  cohort_gradient_richardson = false;
  cohort_gradient_richardson_depth = 4;

  // integration rules for environment splines
  plant_assimilation_adaptive = true;

  plant_assimilation_over_distribution = false;
  plant_assimilation_tol = 1e-6;
  plant_assimilation_iterations = 1000;
  plant_assimilation_rule = 21;

  // used to optimise cohort introduction schedules
  build_schedule_nsteps = 20;
  build_schedule_eps = 1e-3;
  build_schedule_verbose = false;

  // used to find equilibrium seed rain
  equilibrium_nsteps = 20;
  equilibrium_eps = 1e-5;
  equilibrium_large_birth_rate_change = 10;
  equilibrium_verbose = true;
  equilibrium_solver_name = "iteration";
  equilibrium_extinct_birth_rate = 1e-3;
  equilibrium_nattempts = 5;
  equilibrium_solver_logN = true;
  equilibrium_solver_try_keep = true;

  // ODE configuration
  ode_tol_rel = 1e-6;
  ode_tol_abs = 1e-6;
  ode_a_y = 1.0;
  ode_a_dydt = 0.0;
  ode_step_size_min = 1e-6;
  ode_step_size_max = 1e-1;
  ode_step_size_initial = 1e-6;
}

void Control::initialize() {
  if (!plant_assimilation_adaptive) {
    plant_assimilation_iterations = 1;
  }
  integrator =
      quadrature::QAG(plant_assimilation_rule, plant_assimilation_iterations,
                      plant_assimilation_tol, plant_assimilation_tol);
}

} // namespace plant
