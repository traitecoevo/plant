// -*-c++-*-
#ifndef TREE2_CONTROL_H_
#define TREE2_CONTROL_H_

#include <tree2/qag.h>

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
namespace tree2 {

struct Control {
  Control();
  void initialize();

  bool   plant_assimilation_adaptive;

  bool   plant_assimilation_over_distribution;
  double plant_assimilation_tol;
  size_t plant_assimilation_iterations;
  size_t plant_assimilation_rule;
  bool   plant_assimilation_reuse_intervals;

  bool   plant_assimilation_approximate_use;
  double plant_assimilation_approximate_tol;
  int    plant_assimilation_approximate_nbase;
  int    plant_assimilation_approximate_max_depth;
  bool   plant_assimilation_approximate_akima;
  bool   plant_assimilation_approximate_linear;
  bool   plant_assimilation_approximate_rescale_usually;

  double plant_seed_tol;
  int    plant_seed_iterations;

  double cohort_gradient_eps;
  int    cohort_gradient_direction;
  bool   cohort_gradient_richardson;
  size_t cohort_gradient_richardson_depth;

  double environment_light_tol;
  int    environment_light_nbase;
  int    environment_light_max_depth;
  bool   environment_light_akima;
  bool   environment_light_linear;
  bool   environment_light_rescale_usually;
  bool   environment_light_skip;

  double ode_step_size_min;
  double ode_step_size_max;
  double ode_tol_rel;
  double ode_tol_abs;
  double ode_a_y;
  double ode_a_dydt;

  int    schedule_nsteps;
  double schedule_eps;
  bool   schedule_progress;
  bool   schedule_verbose;
  double schedule_default_patch_survival;
  double schedule_default_multipler;
  double schedule_default_min_step_size;
  double schedule_default_max_step_size;

  int    equilibrium_nsteps;
  double equilibrium_eps;
  double equilibrium_large_seed_rain_change;
  bool   equilibrium_progress;
  bool   equilibrium_verbose;
  int    equilibrium_solver;
  double equilibrium_extinct_seed_rain;
  double equilibrium_runsteady_tol;
  double equilibrium_inviable_test_eps;
  int    equilibrium_nattempts;
  bool   equilibrium_solver_logN;
  bool   equilibrium_solver_try_keep;

  // Things derived from this:
  quadrature::QAG integrator;
};

}

#endif
