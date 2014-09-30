// -*-c++-*-
#ifndef TREE_CONTROL_H_
#define TREE_CONTROL_H_

#include <Rcpp.h>

#include "lookup.h"
#include "util.h"
#include "ode_control.h"

namespace model {

class Control : public util::Lookup {
public:
  typedef util::PtrWrapper<Control> ptr;
  Control();
  Control(Rcpp::List x);

  // * R interface
  ode::OdeControl r_ode_control() const;

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

  ode::OdeControl ode_control;

private:
  void do_build_lookup();
  void reset();
  void set_parameters_post_hook();
  ode::OdeControl make_ode_control();

  double _plant_assimilation_adaptive;

  double _plant_assimilation_over_distribution;
  double _plant_assimilation_iterations;
  double _plant_assimilation_rule;
  double _plant_assimilation_reuse_intervals;

  double _plant_assimilation_approximate_use;
  double _plant_assimilation_approximate_nbase;
  double _plant_assimilation_approximate_max_depth;
  double _plant_assimilation_approximate_akima;
  double _plant_assimilation_approximate_linear;
  double _plant_assimilation_approximate_rescale_usually;

  double _plant_seed_iterations;

  double _cohort_gradient_direction;
  double _cohort_gradient_richardson;
  double _cohort_gradient_richardson_depth;

  double _environment_light_nbase;
  double _environment_light_max_depth;
  double _environment_light_akima;
  double _environment_light_linear;
  double _environment_light_rescale_usually;
  double _environment_light_skip;

  double _schedule_nsteps;
  double _schedule_progress;
  double _schedule_verbose;

  double _equilibrium_nsteps;
  double _equilibrium_progress;
  double _equilibrium_verbose;
  double _equilibrium_solver;
};

}

RCPP_EXPOSED_CLASS_NODECL(model::Control)

#endif
