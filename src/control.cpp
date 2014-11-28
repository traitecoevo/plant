#include <tree/control.h>

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

  schedule_nsteps   = 20;
  schedule_eps      = 1e-3;
  schedule_progress = false;
  schedule_verbose  = false;
  // These are designed to agree with Daniel's implementation of the
  // model.
  schedule_default_patch_survival = 6.25302620663814e-05;
  schedule_default_multipler     = 0.2;
  schedule_default_min_step_size = 1e-5;
  schedule_default_max_step_size = 2.0;

  equilibrium_nsteps   = 20;
  equilibrium_eps      = 1e-5;
  equilibrium_large_seed_rain_change = 10;
  equilibrium_progress = false;
  equilibrium_verbose  = true;
  equilibrium_solver   = 1; // TODO: This is a hack
  equilibrium_extinct_seed_rain = 1e-3;
  equilibrium_runsteady_tol = 1e-2;
  equilibrium_inviable_test_eps = 1e-2;
  equilibrium_nattempts = 5;
  equilibrium_solver_logN = true;
  equilibrium_solver_try_keep = true;

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

  _schedule_nsteps =
    static_cast<double>(schedule_nsteps);
  _schedule_progress =
    static_cast<double>(schedule_progress);
  _schedule_verbose =
    static_cast<double>(schedule_verbose);

  _equilibrium_nsteps =
    static_cast<double>(equilibrium_nsteps);
  _equilibrium_progress =
    static_cast<double>(equilibrium_progress);
  _equilibrium_verbose =
    static_cast<double>(equilibrium_verbose);
  _equilibrium_solver =
    static_cast<double>(equilibrium_solver);
  _equilibrium_nattempts =
    static_cast<double>(equilibrium_nattempts);
  _equilibrium_solver_logN =
    static_cast<double>(equilibrium_solver_logN);
  _equilibrium_solver_try_keep =
    static_cast<double>(equilibrium_solver_try_keep);

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

  lookup_table["schedule_nsteps"] =
    &_schedule_nsteps;
  lookup_table["schedule_eps"] =
    &schedule_eps;
  lookup_table["schedule_progress"] =
    &_schedule_progress;
  lookup_table["schedule_verbose"] =
    &_schedule_verbose;
  lookup_table["schedule_default_patch_survival"] =
    &schedule_default_patch_survival;
  lookup_table["schedule_default_multipler"] =
    &schedule_default_multipler;
  lookup_table["schedule_default_min_step_size"] =
    &schedule_default_min_step_size;
  lookup_table["schedule_default_max_step_size"] =
    &schedule_default_max_step_size;

  lookup_table["equilibrium_nsteps"] =
    &_equilibrium_nsteps;
  lookup_table["equilibrium_eps"] =
    &equilibrium_eps;
  lookup_table["equilibrium_large_seed_rain_change"] =
    &equilibrium_large_seed_rain_change;
  lookup_table["equilibrium_progress"] =
    &_equilibrium_progress;
  lookup_table["equilibrium_verbose"] =
    &_equilibrium_verbose;
  lookup_table["equilibrium_solver"] =
    &_equilibrium_solver;
  lookup_table["equilibrium_extinct_seed_rain"] =
    &equilibrium_extinct_seed_rain;
  lookup_table["equilibrium_runsteady_tol"] =
    &equilibrium_runsteady_tol;
  lookup_table["equilibrium_inviable_test_eps"] =
    &equilibrium_inviable_test_eps;
  lookup_table["equilibrium_nattempts"] =
    &_equilibrium_nattempts;
  lookup_table["equilibrium_solver_logN"] =
    &_equilibrium_solver_logN;
  lookup_table["equilibrium_solver_try_keep"] =
    &_equilibrium_solver_try_keep;
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

  schedule_nsteps =
    static_cast<int>(_schedule_nsteps);
  schedule_progress =
    static_cast<int>(_schedule_progress);
  schedule_verbose =
    static_cast<int>(_schedule_verbose);

  equilibrium_nsteps =
    static_cast<int>(_equilibrium_nsteps);
  equilibrium_progress =
    static_cast<bool>(_equilibrium_progress);
  equilibrium_verbose =
    static_cast<bool>(_equilibrium_verbose);
  equilibrium_solver =
    static_cast<int>(_equilibrium_solver);
  equilibrium_nattempts =
    static_cast<int>(_equilibrium_nattempts);
  equilibrium_solver_logN =
    static_cast<int>(_equilibrium_solver_logN);
  equilibrium_solver_try_keep =
    static_cast<int>(_equilibrium_solver_try_keep);

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
