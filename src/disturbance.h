// -*-c++-*-
#ifndef TREE_DISTURBANCE_H_
#define TREE_DISTURBANCE_H_

#include "util.h"
#include "lookup.h"

namespace model {

// Currently hard coded for Weibull only.  Alternatives are
// exponential waiting times, something with an ODE, or something that
// depends on the state.

class Disturbance : public util::Lookup {
public:
  Disturbance();
  double survival_probability(double time_start, double time) const;
private:
  double survival0(double time) const;

  // Lookup related things:
  void do_build_lookup();
  void reset();
  void set_parameters_post_hook();

  double shape;
  double mean_disturbance_interval;
  double scale;
  double p0;
};

}

RCPP_EXPOSED_CLASS(model::Disturbance)

#endif
