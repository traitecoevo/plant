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
  Disturbance(double mean_disturbance_interval);
  double survival_probability(double time_start, double time) const;
private:
  double survival0(double time) const;

  const double shape;
  const double scale;
  const double p0;
};

}
#endif
