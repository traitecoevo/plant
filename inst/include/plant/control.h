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
namespace plant {
struct Control {
  Control();

  size_t function_integration_rule;

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

  size_t schedule_nsteps;
  double schedule_eps;
  bool   schedule_verbose;

  size_t equilibrium_nsteps;
  double equilibrium_eps;
  double equilibrium_large_birth_rate_change;
  bool   equilibrium_verbose;
  std::string equilibrium_solver_name;
  double equilibrium_extinct_birth_rate;
  size_t equilibrium_nattempts;
  bool   equilibrium_solver_logN;
  bool   equilibrium_solver_try_keep;

  bool   save_RK45_cache;
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
