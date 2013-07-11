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
  Disturbance(double mean_interval_);
  double survival_probability(double time_start, double time) const;
  double density(double time) const;
  double r_mean_interval() const;
private:
  double survival0(double time) const;

  double shape;
  double mean_interval;
  double scale;
  double p0;
};

}

RCPP_EXPOSED_CLASS(model::Disturbance)

#endif
