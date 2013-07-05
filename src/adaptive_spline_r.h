// -*-c++-*-
#ifndef TREE_ADAPTIVE_SPLINE_R_H_
#define TREE_ADAPTIVE_SPLINE_R_H_

#include <Rcpp.h>

#include "adaptive_spline.h"
#include "functor.h"

namespace spline {

// This is used for testing a few things, and belongs in the util
// namespace, I think.
class RFunctionWrapper {
public:
  RFunctionWrapper(SEXP fun_, SEXP env_) : fun(fun_), env(env_) {}
  double target(double x) {
    return REAL(Rf_eval(Rf_lang2(fun, Rf_ScalarReal(x)), env))[0];
  }
private:
  SEXP fun, env;
};

namespace test {

Spline test_adaptive_spline(SEXP fun, SEXP env,
			    double a, double b,
			    bool akima_spline);
}


}

#endif
