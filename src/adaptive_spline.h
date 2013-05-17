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
  AdaptiveSpline(util::DFunctor *target);
  void set_control(double atol_, double rtol_,
		   int nbase_, int max_depth_);
  double eval_target(double x) const;
  Spline construct_spline(double a, double b);
private:
  void compute_spline();
  bool refine();

  bool check_err(double y_true, double y_pred) const;
  void check_bounds(double a, double b);

  util::DFunctor *target;

  // Control parameters:
  double atol, rtol, dx, dxmin;
  int nbase, max_depth;

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
