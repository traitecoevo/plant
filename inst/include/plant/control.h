// -*-c++-*-
#ifndef PLANT_PLANT_CONTROL_H_
#define PLANT_PLANT_CONTROL_H_

#include <plant/qag.h>
#include <plant/ode_control.h>
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
//
// TODO: Eventually I need to make sure that the numbers here are
// reasonable, and probably shepherd the translation from int to
// size_t.
namespace plant {

struct Control {
  Control();
  void initialize();

  // patch settings
  double patch_area;
  double patch_max_lifetime; // designed to agree with Daniel's implementation

  // integration rules for net_reproduction
  double plant_seed_tol;
  size_t plant_seed_iterations;

  // integration rules for cohorts
  double cohort_gradient_eps;
  int    cohort_gradient_direction;
  bool   cohort_gradient_richardson;
  size_t cohort_gradient_richardson_depth;

  // integration rules for environment splines (I think..)
  bool   plant_assimilation_adaptive;

  bool   plant_assimilation_over_distribution;
  double plant_assimilation_tol;
  size_t plant_assimilation_iterations;
  size_t plant_assimilation_rule;

  // used to optimise cohort introduction schedules
  size_t build_schedule_nsteps;
  double build_schedule_eps;
  bool   build_schedule_verbose;

  // used to find equilibrium seed rain
  size_t equilibrium_nsteps;
  double equilibrium_eps;
  double equilibrium_large_birth_rate_change;
  bool   equilibrium_verbose;
  std::string equilibrium_solver_name;
  double equilibrium_extinct_birth_rate;
  size_t equilibrium_nattempts;
  bool   equilibrium_solver_logN;
  bool   equilibrium_solver_try_keep;

  // ODE configuration
  double ode_tol_rel;
  double ode_tol_abs;
  double ode_a_y;
  double ode_a_dydt;
  double ode_step_size_min;
  double ode_step_size_max;
  double ode_step_size_initial;

  // Things derived from this:
  quadrature::QAG integrator;
};

inline ode::OdeControl make_ode_control(const Control& control) {
  return ode::OdeControl(control.ode_tol_abs,
                         control.ode_tol_rel,
                         control.ode_a_y,
                         control.ode_a_dydt,
                         control.ode_step_size_min,
                         control.ode_step_size_max,
                         control.ode_step_size_initial);
}

}

#endif
