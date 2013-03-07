// -*-c++-*-
#include <Rcpp.h>

#include "AdaptiveSpline.h"

class RAdaptiveSpline : public AdaptiveSpline {
public:
  RAdaptiveSpline(SEXP fun, SEXP env, double a, double b);
  double target(double x);
private:
  SEXP fun, env;
};
