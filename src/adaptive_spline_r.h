// -*-c++-*-
#ifndef TREE_ADAPTIVE_SPLINE_R_H_
#define TREE_ADAPTIVE_SPLINE_R_H_

#include <Rcpp.h>

#include "adaptive_spline.h"

namespace spline {

class AdaptiveSplineR : public AdaptiveSpline {
public:
  AdaptiveSplineR(SEXP fun, SEXP env, double a, double b);
  double target(double x);
private:
  SEXP fun, env;
};

}

#endif
