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

  // Number of collocation nodes for the multirate (method="mri") fast sub-cycle.
  // The soil sub-cycle needs the per-layer root uptake at each micro-step; that
  // uptake is a density-weighted integral of per-cohort consumption over the
  // size distribution. 0 (the default) evaluates it over all N cohorts (exact);
  // m > 0 quadratures it at m of the N frozen cohorts instead, so the sub-cycle
  // costs m physiology solves, not N. Converges ~O(m^-2), but the accuracy at a
  // given m depends strongly on the distribution: young stands reach the sub-1%
  // range by m≈20, mature (skewed) stands need m close to N for the same
  // accuracy under the current even-index node placement (a smarter,
  // importance-weighted placement is the open item -- plant#53 §6/§4.4). Only
  // consulted on the method="mri" path; ignored otherwise.
  size_t n_collocation_nodes;

  // Multirate (method="mri") fast-block inner stepper. false (default) sub-cycles
  // the soil block with the adaptive black-box RK; true uses the exact-flow split
  // (R1 analytic drainage recession + ROS34PW2 on the gentle remainder), which
  // removes the drainage stiffness so the sub-cycle takes far fewer micro steps
  // (Lever 1). Only consulted on the method="mri" path.
  bool mri_use_split;

  // ODE integration method for the SCM resident solver. One of "rkck" (the
  // default adaptive Cash-Karp explicit RK), "rodas" (the stiff Rosenbrock
  // stepper), or "mri" (the multirate MRI-GARK stepper: a fixed macro grid that
  // sub-cycles the fast soil column, for TF24). Empty is treated as "rkck", so
  // default behaviour is unchanged. Selected once when the SCM builds its
  // Solver; every other integration path is untouched.
  std::string ode_method;

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
