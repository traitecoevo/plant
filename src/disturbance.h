// -*-c++-*-
#ifndef TREE_DISTURBANCE_H_
#define TREE_DISTURBANCE_H_

#include "util.h"

namespace model {

// Currently hard coded for Weibull only.  Alternatives are
// exponential waiting times, something with an ODE, or something that
// depends on the state.

class Disturbance {
public:
  typedef util::PtrWrapper<Disturbance> ptr;
  // TODO: 2nd contructor would be better if this was util::Lookup?
  Disturbance();
  Disturbance(double mean_disturbance_interval);
  double survival_probability(double time_start, double time) const;
private:
  double survival0(double time) const;
  void set_mean_disturbance_interval(double x);

  const double shape;
  double scale;
  double p0;
};

}

RCPP_EXPOSED_CLASS(model::Disturbance)

#endif
