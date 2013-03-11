#include "adaptive_spline_r.h"

namespace spline {

AdaptiveSplineR::AdaptiveSplineR(SEXP fun_, SEXP env_, 
				 double a, double b) :
  AdaptiveSpline(), fun(fun_), env(env_) {
  set_bounds(a, b);
  set_target(helper_spline<AdaptiveSplineR>, this);
  construct_spline();
}

double AdaptiveSplineR::target(double x) {
  return REAL(Rf_eval(Rf_lang2(fun, Rf_ScalarReal(x)), env))[0];
}

}
