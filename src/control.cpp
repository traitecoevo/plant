#include <plant/control.h>

namespace plant {

Control::Control() {

  // Number of points used when numerically intergrating a function
  // using Gauss-Kronrod quadrature. Rules defined in qk_rules.cpp
  function_integration_rule = 21;

  offspring_production_tol= 1e-8;
  offspring_production_iterations = 1000;

  node_gradient_eps = 1e-6;
  node_gradient_direction = 1;
  node_gradient_richardson = false;
  node_gradient_richardson_depth = 4;

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

  save_RK45_cache = false;
}

}
