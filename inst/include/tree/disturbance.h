// -*-c++-*-
#ifndef TREE_TREE_DISTURBANCE_H_
#define TREE_TREE_DISTURBANCE_H_

namespace tree {

// Currently hard coded for Weibull only.  Alternatives are
// exponential waiting times, something with an ODE, or something that
// depends on the state.

class Disturbance {
public:
  Disturbance(double mean_interval_);
  double density(double time) const;
  double r_mean_interval() const;
  double pr_survival(double time) const;
  double pr_survival_conditional(double time, double time_start) const;
  double cdf(double pr) const;

private:
  double shape;
  double mean_interval;
  double scale;
  double p0;
};


}

#endif
