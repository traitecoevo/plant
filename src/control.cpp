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

  // Crown shading model (see Control header). Empty = each strategy's own
  // default (FF16 -> deep-crown, TF24 -> mean-light), so default behaviour is
  // unchanged for both.
  shading_model = "";

  // PPA canopy layer thickness in optical-depth units (see Control header).
  ppa_layer_optical_depth = 0.5;
  // PPA layer-boundary smoothing fraction (see Control header).
  ppa_layer_smoothing = 0.3;

  offspring_production_tol= 1e-8;
  offspring_production_iterations = 1000;

  node_gradient_eps = 1e-6;
  node_gradient_direction = -1;
  node_gradient_richardson = false;
  node_gradient_richardson_depth = 4;
  node_gradient_exact_ad = false;

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
  schedule_eps      = 2e-2;
  schedule_verbose  = false;

  save_RK45_cache = false;

  GSS_tol_abs = 1e-3;
  vulnerability_curve_ncontrol = 1e2;
  ci_abs_tol = 1e-3;
  ci_niter = 1e3;
}

}
