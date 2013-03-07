#ifndef RADAPTIVE_SPLINE_H
#define RADAPTIVE_SPLINE_H

#include "RAdaptiveSpline.h"

RAdaptiveSpline::RAdaptiveSpline(SEXP fun_, SEXP env_, 
				 double a, double b) :
  AdaptiveSpline(), fun(fun_), env(env_) {
  set_bounds(a, b);
  set_target(helper_spline<RAdaptiveSpline>, this);
  construct_spline();
}

double RAdaptiveSpline::target(double x) {
  return REAL(Rf_eval(Rf_lang2(fun, Rf_ScalarReal(x)), env))[0];
}

#endif
