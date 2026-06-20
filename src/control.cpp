#include <plant/control.h>

namespace plant {

Control::Control() {
  // These defaults are the pragmatic "fast-ish" settings used for essentially
  // all of plant's runs. They previously lived in the R helpers fast_control()
  // / scm_base_control(), which layered them on top of a tighter-tolerance
  // constructor; those helpers have been removed and the values folded in here
  // so Control() is the single source of truth. Tighten tolerances explicitly
  // (e.g. ode_tol_rel/abs, schedule_eps) if you need a high-accuracy run.

  // Number of points used when numerically intergrating a function
  // using Gauss-Kronrod quadrature. Rules defined in qk_rules.cpp
  function_integration_rule = 21;

  offspring_production_tol= 1e-8;
  offspring_production_iterations = 1000;

  node_gradient_eps = 1e-6;
  node_gradient_direction = -1;
  node_gradient_richardson = false;
  node_gradient_richardson_depth = 4;

  ode_step_size_initial = 1e-6;
  ode_step_size_min = 1e-6;
  ode_step_size_max = 5;
  ode_tol_rel       = 1e-4;
  ode_tol_abs       = 1e-4;
  ode_a_y           = 1.0;
  ode_a_dydt        = 0.0;

  // 0 = adaptive RKCK (default); > 0 selects fixed-step forward Euler with this
  // spacing in years (see control.h).
  fixed_time_step   = 0.0;

  schedule_nsteps   = 20;
  schedule_eps      = 5e-3;
  schedule_verbose  = false;

  save_RK45_cache = false;

  GSS_tol_abs = 1e-3;
  vulnerability_curve_ncontrol = 1e2;
  ci_abs_tol = 1e-3;
  ci_niter = 1e3;
}

}
