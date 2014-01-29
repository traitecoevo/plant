#include "control.h"

namespace model {

Control::Control() {
  reset();
}

Control::Control(Rcpp::List x) {
  reset();
  set_parameters(x);
}

// * R interface
ode::OdeControl Control::r_ode_control() const {
  return ode_control;
}

// * Private methods
void Control::reset() {
  plant_assimilation_adaptive = true;

  plant_assimilation_over_distribution = false;
  plant_assimilation_tol = 1e-6;
  plant_assimilation_iterations = 1000;
  plant_assimilation_rule = 21;
  plant_assimilation_reuse_intervals = true;

  plant_assimilation_approximate_use = false;
  plant_assimilation_approximate_tol = 1e-6;
  plant_assimilation_approximate_nbase = 17;
  plant_assimilation_approximate_max_depth = 16;
  plant_assimilation_approximate_akima = false;
  plant_assimilation_approximate_linear = false;
  plant_assimilation_approximate_rescale_usually = false;

  plant_seed_tol = 1e-6;
  plant_seed_iterations = 1000;

  cohort_gradient_eps = 1e-6;
  cohort_gradient_direction = 1;
  cohort_gradient_richardson = false;
  cohort_gradient_richardson_depth = 4;

  environment_light_tol = 1e-6;
  environment_light_nbase = 17;
  environment_light_max_depth = 16;
  environment_light_akima = false;
  environment_light_linear = false;
  environment_light_rescale_usually = false;
  environment_light_skip = false;

  ode_step_size_min = 1e-6;
  ode_step_size_max = 1e-1;
  ode_tol_rel       = 1e-6;
  ode_tol_abs       = 1e-6;
  ode_a_y           = 1.0;
  ode_a_dydt        = 0.0;
  // TODO: Also no_steps_max?

  // Then set the values for the lookup table, based on these (this is
  // basically the inverse of set_parameters_post_hook())
  _plant_assimilation_adaptive =
    static_cast<double>(plant_assimilation_adaptive);

  _plant_assimilation_over_distribution =
    static_cast<double>(plant_assimilation_over_distribution);
  _plant_assimilation_iterations =
    static_cast<double>(plant_assimilation_iterations);
  _plant_assimilation_rule =
    static_cast<double>(plant_assimilation_rule);
  _plant_assimilation_reuse_intervals =
    static_cast<double>(plant_assimilation_reuse_intervals);

  _plant_assimilation_approximate_use =
    static_cast<bool>(plant_assimilation_approximate_use);
  _plant_assimilation_approximate_nbase =
    static_cast<double>(plant_assimilation_approximate_nbase);
  _plant_assimilation_approximate_max_depth =
    static_cast<double>(plant_assimilation_approximate_max_depth);
  _plant_assimilation_approximate_akima =
    static_cast<bool>(plant_assimilation_approximate_akima);
  _plant_assimilation_approximate_linear =
    static_cast<bool>(plant_assimilation_approximate_linear);
  _plant_assimilation_approximate_rescale_usually =
    static_cast<bool>(plant_assimilation_approximate_rescale_usually);

  _plant_seed_iterations =
    static_cast<double>(plant_seed_iterations);

  _cohort_gradient_direction =
    static_cast<double>(cohort_gradient_direction);
  _cohort_gradient_richardson =
    static_cast<double>(cohort_gradient_richardson);
  _cohort_gradient_richardson_depth =
    static_cast<double>(cohort_gradient_richardson_depth);

  _environment_light_nbase =
    static_cast<double>(environment_light_nbase);
  _environment_light_max_depth =
    static_cast<double>(environment_light_max_depth);
  _environment_light_akima =
    static_cast<bool>(environment_light_akima);
  _environment_light_linear =
    static_cast<bool>(environment_light_linear);
  _environment_light_rescale_usually =
    static_cast<bool>(environment_light_rescale_usually);
  _environment_light_skip =
    static_cast<bool>(environment_light_skip);

  // Like set_parameters_post_hook(), rebuild the ODE control, too.
  ode_control = make_ode_control();
}

void Control::do_build_lookup() {
  lookup_table["plant_assimilation_adaptive"] =
    &_plant_assimilation_adaptive;

  lookup_table["plant_assimilation_over_distribution"] =
    &_plant_assimilation_over_distribution;
  lookup_table["plant_assimilation_tol"] =
    &plant_assimilation_tol;
  lookup_table["plant_assimilation_iterations"] =
    &_plant_assimilation_iterations;
  lookup_table["plant_assimilation_rule"] =
    &_plant_assimilation_rule;
  lookup_table["plant_assimilation_reuse_intervals"] =
    &_plant_assimilation_reuse_intervals;

  lookup_table["plant_assimilation_approximate_use"] =
    &_plant_assimilation_approximate_use;
  lookup_table["plant_assimilation_approximate_tol"] =
    &plant_assimilation_approximate_tol;
  lookup_table["plant_assimilation_approximate_nbase"] =
    &_plant_assimilation_approximate_nbase;
  lookup_table["plant_assimilation_approximate_max_depth"] =
    &_plant_assimilation_approximate_max_depth;
  lookup_table["plant_assimilation_approximate_akima"] =
    &_plant_assimilation_approximate_akima;
  lookup_table["plant_assimilation_approximate_linear"] =
    &_plant_assimilation_approximate_linear;
  lookup_table["plant_assimilation_approximate_rescale_usually"] =
    &_plant_assimilation_approximate_rescale_usually;

  lookup_table["plant_seed_tol"] =
    &plant_seed_tol;
  lookup_table["plant_seed_iterations"] =
    &_plant_seed_iterations;

  lookup_table["cohort_gradient_eps"] = 
    &cohort_gradient_eps;
  lookup_table["cohort_gradient_direction"] =
    &_cohort_gradient_direction;
  lookup_table["cohort_gradient_richardson"] =
    &_cohort_gradient_richardson;
  lookup_table["cohort_gradient_richardson_depth"] =
    &_cohort_gradient_richardson_depth;

  lookup_table["environment_light_tol"] =
    &environment_light_tol;
  lookup_table["environment_light_nbase"] =
    &_environment_light_nbase;
  lookup_table["environment_light_max_depth"] =
    &_environment_light_max_depth;
  lookup_table["environment_light_akima"] =
    &_environment_light_akima;
  lookup_table["environment_light_linear"] =
    &_environment_light_linear;
  lookup_table["environment_light_rescale_usually"] =
    &_environment_light_rescale_usually;
  lookup_table["environment_light_skip"] =
    &_environment_light_skip;

  lookup_table["ode_step_size_min"] =
    &ode_step_size_min;
  lookup_table["ode_step_size_max"] =
    &ode_step_size_max;
  lookup_table["ode_tol_rel"] =
    &ode_tol_rel;
  lookup_table["ode_tol_abs"] =
    &ode_tol_abs;
  lookup_table["ode_a_y"] =
    &ode_a_y;
  lookup_table["ode_a_dydt"] =
    &ode_a_dydt;
}

void Control::set_parameters_post_hook() {
  plant_assimilation_adaptive =
    _plant_assimilation_adaptive != 0.0;

  plant_assimilation_over_distribution =
    _plant_assimilation_over_distribution != 0.0;
  plant_assimilation_iterations =
    static_cast<size_t>(_plant_assimilation_iterations);
  plant_assimilation_rule =
    static_cast<size_t>(_plant_assimilation_rule);
  plant_assimilation_reuse_intervals =
    _plant_assimilation_reuse_intervals != 0.0;

  plant_assimilation_approximate_use =
    static_cast<bool>(_plant_assimilation_approximate_use);
  plant_assimilation_approximate_nbase =
    static_cast<int>(_plant_assimilation_approximate_nbase);
  plant_assimilation_approximate_max_depth =
    static_cast<int>(_plant_assimilation_approximate_max_depth);
  plant_assimilation_approximate_akima =
    static_cast<bool>(_plant_assimilation_approximate_akima);
  plant_assimilation_approximate_linear =
    static_cast<bool>(_plant_assimilation_approximate_linear);
  plant_assimilation_approximate_rescale_usually =
    static_cast<bool>(_plant_assimilation_approximate_rescale_usually);

  plant_seed_iterations =
    static_cast<int>(_plant_seed_iterations);

  cohort_gradient_direction =
    static_cast<int>(_cohort_gradient_direction);
  cohort_gradient_richardson = _cohort_gradient_richardson != 0.0;
  cohort_gradient_richardson_depth =
    static_cast<size_t>(_cohort_gradient_richardson_depth);

  environment_light_nbase =
    static_cast<int>(_environment_light_nbase);
  environment_light_max_depth =
    static_cast<int>(_environment_light_max_depth);
  environment_light_akima =
    static_cast<bool>(_environment_light_akima);
  environment_light_linear =
    static_cast<bool>(_environment_light_linear);
  environment_light_rescale_usually =
    static_cast<bool>(_environment_light_rescale_usually);
  environment_light_skip =
    static_cast<bool>(_environment_light_skip);

  ode_control = make_ode_control();
}

ode::OdeControl Control::make_ode_control() {
  ode::OdeControl ret(ode_tol_rel,
		      ode_tol_abs,
		      ode_a_y,
		      ode_a_dydt,
		      ode_step_size_min,
		      ode_step_size_max);
  return ret;
}

}
