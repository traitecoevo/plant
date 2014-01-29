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
  void set_time(double t);

  double eval(double u)  const {return current.eval(u);}
  double deriv(double u) const {return current.deriv(u);}
  double min()           const {return current.min();}
  double max()           const {return current.max();}

  Interpolator get_current() const {return current;}

  Interpolator merge(double t,
		     double t0, const Interpolator& s0,
		     double t1, const Interpolator& s1) const;

private:
  void initialise();

  std::vector<double> times;
  std::vector<Interpolator> splines;
  double time;
  double time_max;
  Interpolator current;
};

}

#endif
