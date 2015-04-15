// -*-c++-*-
#ifndef TREE_FAKE_LIGHT_ENVIRONMENT_H_
#define TREE_FAKE_LIGHT_ENVIRONMENT_H_

// The name here will probably change; this is a quick and dirty 2d
// interpolator that is tuned to a specific problem though.
//
// There are two possible solutions; we can either lookup and
// interpolate (that is, keep track of two splines - lower and upper -
// evaluate both splines and interpolate between those points) or we
// can interpolate and lookup (that is, create a new spline from the
// interpolatoion of two splines pair of splines).  I'm going to try
// and create a spline.

#include <Rcpp.h>
#include "interpolator.h"

namespace interpolator {

class FakeLightEnvironment {
public:
  FakeLightEnvironment(std::vector<double> times_,
		       std::vector<Interpolator> splines_);
  FakeLightEnvironment(std::vector<double> times_,
		       Rcpp::List splines_);
  Interpolator operator()(double t);

private:
  void initialise();
  Interpolator merge(double t,
		     double t0, const Interpolator& s0,
		     double t1, const Interpolator& s1) const;

  std::vector<double> times;
  std::vector<Interpolator> splines;
  double time_max;
};

}

RCPP_EXPOSED_CLASS_NODECL(interpolator::FakeLightEnvironment)

#endif
