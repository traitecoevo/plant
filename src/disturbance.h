// -*-c++-*-
#ifndef TREE_DISTURBANCE_H_
#define TREE_DISTURBANCE_H_

namespace model {

// Currently hard coded for Weibull only.  Alternatives are
// exponential waiting times, something with an ODE, or something that
// depends on the state.

class Disturbance {
public:
  Disturbance(double mean_disturbance_interval);
  double survival_probability(double time_start, double time);
private:
  double survival0(double time);

  const double shape;
  const double scale;
  const double p0;
};

}
#endif
