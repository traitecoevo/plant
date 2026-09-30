// -*-c++-*-
#ifndef PLANT_PLANT_CONTROL_H_
#define PLANT_PLANT_CONTROL_H_

#include <plant/qag.h>
#include <odelia/ode_control.hpp>
#include <string>

// The `Control` object holds all the non-biological control
// parameters.  These might get templated against different ways of
// running things as in the same way as `Strategy` and `Parameters`
// but for now assume that they don't.
//
// Control is really hierarchical, but is not actually modelled that
// way yet.  For now, the hierarchy is indicated only by naming
// convention, but this is stored as a flat bunch of things.
//
// Because Control is essentially a dumb set of parameters that has no
// real functionality, we don't export it as a reference class, but
// instead use RcppR6's "list" export ability.
namespace plant {
struct Control {
  Control();

  size_t function_integration_rule;

  // Crown shading model. One of "deep-crown" (integrate photosynthesis over
  // crown depth), "mean-light" (integrate the light over depth, one
  // photosynthesis evaluation of the mean), "crown-centre" (single evaluation at
  // the crown centre), or "ppa" (FF16 only: discrete stepped light layers). The
  // empty default "" means "use the strategy's own default": FF16 -> deep-crown,
  // TF24 -> mean-light. Resolved once in each strategy's prepare_strategy(),
  // so it never costs a string comparison on the hot path.
  std::string shading_model;

  // PPA only: thickness of one discrete canopy layer, in optical-depth units
  // (tau = sum of k * leaf-area-index above a height). The stepped light
  // profile floors tau to integer multiples of this value. The default 0.5
  // corresponds to one unit of leaf area index per layer at the FF16 default
  // light-extinction coefficient k_I = 0.5. Ignored by the other models.
  double ppa_layer_optical_depth;

  // PPA only: smoothing width of each layer boundary, as a fraction (0, 1] of
  // the layer thickness. The stepped profile is flat over the lower part of
  // each layer and ramps smoothly (cubic smoothstep) over the top `fraction` of
  // it, so the profile is C1-continuous and can be integrated by the adaptive
  // ODE solver. -> 0 approaches a hard step (and its numerical instability);
  // = 1 removes the flat region (approaches the smooth deep-crown profile).
  // This one setting distinguishes the two PPA variants: "PPA (hard step)"
  // (= 0, the literal field discretisation; discontinuous, does not run in the
  // adaptive solver) vs "PPA (smoothed)" (> 0, default 0.3; the runnable
  // version). They are the same `ppa` shading model, not separate models.
  double ppa_layer_smoothing;

  double offspring_production_tol;
  size_t offspring_production_iterations;

  double node_gradient_eps;
  int    node_gradient_direction;
  bool   node_gradient_richardson;
  size_t node_gradient_richardson_depth;

  // Which coordinate the size distribution's density is carried in:
  //
  //   "auto"       let the model decide (the default)
  //   "birth_date" density in birth date; the two resource integrals run over
  //                introduction times and the density rate is mortality alone
  //   "height"     density in height; the density rate additionally carries a
  //                compression term
  //
  // The height coordinate's compression term is the total derivative of the
  // growth rate along a cohort's own trajectory, which equals the partial
  // d(growth)/d(height) only when growth is a function of size alone. A model
  // carrying state that feeds back into growth breaks that, so the coordinate
  // is not a free numerical choice for every model -- hence "auto", resolved
  // per model against Strategy::density_in_birth_date_default. An explicit
  // value is an override and always wins, so both coordinates stay available
  // for the comparison tests. See notes/516-flexible-allometry.md.
  std::string node_density_coordinate;

  double ode_step_size_initial;
  double ode_step_size_min;
  double ode_step_size_max;
  double ode_tol_rel;
  double ode_tol_abs;
  double ode_a_y;
  double ode_a_dydt;

  // Fixed-step ODE integration (forward Euler).  Units: years.  When 0 (the
  // default) the SCM integrates residents with the adaptive, error-controlled
  // Cash-Karp RKCK solver.  When > 0 it instead uses plain forward Euler on a
  // uniform grid of this spacing (e.g. 1/365 for a daily step), the way
  // industry-standard DGVMs are run.  Note: forward Euler is incompatible with
  // the mutant-fitness replay path and with save_RK45_cache (the RK sub-step
  // cache has no Euler analogue); the SCM errors clearly if combined.
  double fixed_time_step;

  size_t schedule_nsteps;
  double schedule_eps;
  bool   schedule_verbose;

  bool   save_RK45_cache;

    //TF24 control parameters
  double GSS_tol_abs;
  double vulnerability_curve_ncontrol;
  double ci_abs_tol;
  double ci_niter;
};

// Resolve Control::node_density_coordinate for a model whose own preference is
// `model_prefers_birth_date` (Strategy::density_in_birth_date_default). Only
// "auto" consults the model; an explicit setting wins. An unrecognised value is
// refused rather than silently treated as one of the two coordinates, because
// the coordinates do not agree for a model with carried state -- a typo would
// otherwise read as a real result.
inline bool density_in_birth_date(const Control& control,
                                  bool model_prefers_birth_date) {
  const std::string& coord = control.node_density_coordinate;
  if (coord == "auto" || coord.empty()) {
    return model_prefers_birth_date;
  }
  if (coord == "birth_date") {
    return true;
  }
  if (coord == "height") {
    return false;
  }
  util::stop("Unknown node_density_coordinate \"" + coord +
             "\": expected \"auto\", \"birth_date\" or \"height\"");
}

inline odelia::ode::OdeControl make_ode_control(const Control& control) {
  return odelia::ode::OdeControl(control.ode_tol_abs,
                         control.ode_tol_rel,
                         control.ode_a_y,
                         control.ode_a_dydt,
                         control.ode_step_size_min,
                         control.ode_step_size_max,
                         control.ode_step_size_initial);
}

}

#endif
