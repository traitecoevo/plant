#include "control.h"

namespace model {

Control::Control() {
  reset();
  set_parameters_post_hook();
}

Control::Control(Rcpp::List x) {
  reset();
  set_parameters(x);
}

void Control::reset() {
  _plant_assimilation_over_distribution = static_cast<double>(false);
  plant_assimilation_tol = 1e-6;
  _plant_assimilation_iterations = static_cast<double>(1000);

  cohort_gradient_eps = 1e-6;
  _cohort_gradient_richardson = static_cast<double>(false);
  _cohort_gradient_richardson_depth = static_cast<double>(4);
}

void Control::do_build_lookup() {
  lookup_table["plant_assimilation_over_distribution"] =
    &_plant_assimilation_over_distribution;
  lookup_table["plant_assimilation_tol"] =
    &plant_assimilation_tol;
  lookup_table["plant_assimilation_iterations"] =
    &_plant_assimilation_iterations;

  lookup_table["cohort_gradient_eps"] = 
    &cohort_gradient_eps;
  lookup_table["cohort_gradient_richardson"] =
    &_cohort_gradient_richardson;
  lookup_table["cohort_gradient_richardson_depth"]=
    &_cohort_gradient_richardson_depth;
}

void Control::set_parameters_post_hook() {
  plant_assimilation_over_distribution =
    _plant_assimilation_over_distribution != 0.0;
  plant_assimilation_iterations =
    static_cast<int>(_plant_assimilation_iterations);

  cohort_gradient_richardson = _cohort_gradient_richardson != 0.0;
  cohort_gradient_richardson_depth =
    static_cast<int>(_cohort_gradient_richardson_depth);
}

}
