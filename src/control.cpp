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
  cohort_gradient_eps = 1e-6;
  _cohort_gradient_richardson = 0.0;
  _cohort_gradient_richardson_depth = 4.0;
}

void Control::do_build_lookup() {
  lookup_table["cohort_gradient_eps"] = &cohort_gradient_eps;
  lookup_table["cohort_gradient_richardson"] =
    &_cohort_gradient_richardson;
  lookup_table["cohort_gradient_richardson_depth"]=
    &_cohort_gradient_richardson_depth;
}

void Control::set_parameters_post_hook() {
  cohort_gradient_richardson = _cohort_gradient_richardson != 0.0;
  cohort_gradient_richardson_depth =
    static_cast<int>(_cohort_gradient_richardson_depth);
}

}
