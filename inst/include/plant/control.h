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

  // The absorbing top layer of a continuous patch (#519): with epsilon > 0, a
  // cohort that has grown to the height where its strategy allocates 1 - epsilon
  // of its maximum to reproduction (Strategy::asymptotic_height) is absorbed at
  // the next node introduction into one top cohort per species, pinned at that
  // height, which carries its density, its future seed output and its share of
  // the light field and uptake. Without it a stand that runs on never
  // equilibrates: the oldest cohorts edge up toward their asymptote for as long
  // as they carry density. Zero (the default) absorbs nothing. Birth-date
  // coordinate only, and only for a strategy that declares an asymptotic height.
  double absorb_epsilon;

  double ode_step_size_initial;
  double ode_step_size_min;
  double ode_step_size_max;
  double ode_tol_rel;
  double ode_tol_abs;
  double ode_a_y;
  double ode_a_dydt;

  // Substep TF24's soil where a step is past Cash-Karp's stability limit on it
  // or substepping lets it grow. A step that spans a knot of the rain substeps
  // it whatever this is.
  bool ode_soil_substep;
  // The relative tolerance that chooses the soil's substeps (odelia's
  // OdeControl::substep_tol); the absolute is 1e-4 of it.
  double ode_soil_substep_tol;

  // Multiplies the tolerance on TF24's soil layers, so above 1 loosens it.
  double ode_tol_scale_soil;

  // Fixed-step ODE integration (forward Euler).  Units: years.  When 0 (the
  // default) the SCM integrates residents with the adaptive, error-controlled
  // Cash-Karp RKCK solver.  When > 0 it instead uses plain forward Euler on a
  // uniform grid of this spacing (e.g. 1/365 for a daily step), the way
  // industry-standard DGVMs are run.
  double fixed_time_step;

  size_t schedule_nsteps;
  double schedule_eps;
  bool   schedule_verbose;

  // The smallest profit curvature the interior derivation will divide by. Every
  // collar response is a quotient over it, so a curvature approaching zero
  // returns amplification rather than a derivative -- and returns it finite,
  // which is the shape nothing catches. Below this the water rows are refused
  // and the profit row, which reads no curvature, is emitted anyway.
  double gradient_curvature_floor;

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

// Whether a cohort with no density is frozen in place -- no rate evaluation,
// every rate zero -- rather than integrated. In a continuous patch a cohort
// leaves the integrated set in one way only, by passing the absorbing height
// into the top layer, which takes its density with it; a cohort whose density
// is already exactly zero is frozen too, since there is nothing of it to keep.
// Its per-capita output stops where it was. That is zero for a cohort that
// failed to establish, and nothing to speak of for one that died, but a cohort
// BORN with no density -- at a zero birth rate -- would be frozen at birth and
// its per-capita output read zero where it would have accrued a little, so a
// zero birth rate is refused at introduction (Patch::push_nodes). Off
// unless the top layer is on, so any other run integrates its empty cohorts
// exactly as it always has.
//
// There is deliberately no survival floor: retiring a cohort below one discards
// the biomass old trees still hold rather than conserving it (#519).
inline bool freezes_retired(const Control& c) {
  return c.absorb_epsilon > 0.0;
}

// Refuses a tolerance the stepper cannot use: one that is not positive.
void validate_ode_tolerances(const Control& control);

// odelia's step-size rule without its stabilisation term: a stand's step is
// set by accuracy and by the introduction schedule rather than by stability,
// and the term only adds steps there (12-20% of FF16's run time).
inline odelia::ode::OdeControl make_ode_control(const Control& control) {
  odelia::ode::OdeControl ret(control.ode_tol_abs,
                              control.ode_tol_rel,
                              control.ode_a_y,
                              control.ode_a_dydt,
                              control.ode_step_size_min,
                              control.ode_step_size_max,
                              control.ode_step_size_initial);
  ret.set_stabilisation(0.0);
  ret.substep_subsystem = control.ode_soil_substep;
  ret.substep_tol = control.ode_soil_substep_tol;
  return ret;
}

}

#endif
