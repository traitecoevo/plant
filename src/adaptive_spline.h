// -*-c++-*-
#ifndef TREE_ADAPTIVE_SPLINE_H_
#define TREE_ADAPTIVE_SPLINE_H_

#include <list>
#include <Rcpp.h>

#include "spline.h"
#include "functor.h"

namespace spline {

class AdaptiveSpline : public Spline {
public:
  AdaptiveSpline();
  void set_bounds(double a_, double b_);
  void set_target(util::DFunctor *target_);
  void set_control(double atol_, double rtol_, 
		   int nbase_, int max_depth_);
  bool construct_spline();

private:
  bool check_err(double y_true, double y_pred) const;
  void compute_spline();
  bool refine();
  void reset();

  void init(); // disable from base.

  void check_bounds();

  double a, b, atol, rtol, dx, dxmin;
  int nbase, max_depth;

  util::DFunctor *target;

  // In contrast to the Spline's x and y, which are vectors, these are
  // lists so that we can easily add points in the middle of them.
  // 'yr' and 'yp' are y "real" and y "predicted", respectively.
  std::list<double> xx, yy;
  std::list<bool> zz;
};

}

RCPP_EXPOSED_CLASS(spline::AdaptiveSpline)

#endif
