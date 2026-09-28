#include <plant/control.h>
#include <phylloptim/leaf_model.hpp>

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
  node_density_in_birth_date = false;

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


  // Measured rather than chosen: over 5625 solved leaf states spanning stem_c
  // 0.4 to 2.68, stem_b 2.5 to 6.0, beta2 0.5 to 3.0, radiation 15 to 2000 and
  // soil potentials 0.1 to 5.5, every one of the 1351 interior points had a
  // strictly negative curvature and the smallest magnitude was 0.0623. This sits
  // sixty times below that, so it separates a divergence from the range the
  // model occupies rather than narrowing what answers.
  //
  // ⚠️ The MAGNITUDE claim above is what this number rests on, and it holds. The
  // SIGN claim -- "every one strictly negative" -- does not: the century stand in
  // scripts/profile-stand-gradient.R refuses on an interior point whose curvature
  // is +34.414226, and no floor admits that point. The sweep was over static leaf
  // states in the box named above; a stand integrated for 105 years reaches
  // operating points outside it. So the guard's sign limb is live in practice and
  // is not a defensive branch -- which is why it says something different from the
  // floor limb when it fires.
  gradient_curvature_floor = 1e-3;

  // Bracket tolerance of the collar-potential search. It has only to land inside
  // the basin of the Newton polish, which sets the operating point returned.
  GSS_tol_abs = 1e-1;
  // The leaf's own choice, read rather than restated: this was 1e2 while
  // phylloptim's default was sixteen times finer, and since every stand builds its
  // leaf through here, that number was the one the whole reverse sweep ran on.
  vulnerability_curve_ncontrol = phylloptim::Leaf::ncontrol_default;
  ci_abs_tol = 1e-3;
  ci_niter = 1e3;
}

}
