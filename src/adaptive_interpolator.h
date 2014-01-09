// -*-c++-*-
#ifndef TREE_ADAPTIVE_INTERPOLATOR_H_
#define TREE_ADAPTIVE_INTERPOLATOR_H_

#include <list>
#include <Rcpp.h>

#include "interpolator.h"
#include "functor.h"

namespace interpolator {

class AdaptiveInterpolator {
public:
  AdaptiveInterpolator(double atol_, double rtol_,
		       int nbase_, int max_depth_,
		       bool akima, bool linear);
  double eval_target(double x) const;
  Interpolator construct(util::DFunctor *target_, double a, double b);
private:
  void compute();
  bool refine();

  bool check_err(double y_true, double y_pred) const;
  static void check_bounds(double a, double b);

  util::DFunctor *target;

  // Control parameters:
  double atol, rtol;
  int nbase, max_depth;
  double dx, dxmin;

  // In contrast to the Interpolator's x and y, which are vectors, these are
  // lists so that we can easily add points in the middle of them.
  std::list<double> xx, yy;
  std::list<bool> zz;

  // Temporary interpolator object that we build.
  Interpolator interpolator;
};

namespace test {
Interpolator test_adaptive_interpolator(SEXP fun, SEXP env,
					double a, double b,
					bool akima, bool linear);
}

}

RCPP_EXPOSED_CLASS(interpolator::AdaptiveInterpolator)

#endif
