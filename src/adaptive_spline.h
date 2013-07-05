// -*-c++-*-
#ifndef TREE_ADAPTIVE_SPLINE_H_
#define TREE_ADAPTIVE_SPLINE_H_

#include <list>
#include <Rcpp.h>

#include "spline.h"
#include "functor.h"

namespace spline {

class AdaptiveSpline {
public:
  AdaptiveSpline(double atol_, double rtol_,
		 int nbase_, int max_depth_,
		 bool akima_spline);
  double eval_target(double x) const;
  Spline construct_spline(util::DFunctor *target_, double a, double b);
private:
  void compute_spline();
  bool refine();

  bool check_err(double y_true, double y_pred) const;
  static void check_bounds(double a, double b);

  util::DFunctor *target;

  // Control parameters:
  double atol, rtol;
  int nbase, max_depth;
  double dx, dxmin;

  // In contrast to the Spline's x and y, which are vectors, these are
  // lists so that we can easily add points in the middle of them.
  std::list<double> xx, yy;
  std::list<bool> zz;

  // Temporary spline object that we build.
  Spline spline;
};

}

RCPP_EXPOSED_CLASS(spline::AdaptiveSpline)

#endif
