#include "disturbance.h"

#include <Rcpp.h>

namespace model {

// This should be rewritten to more closely use R's actual Weibull
// functions.

Disturbance::Disturbance()
  : shape(2.0) {
  set_mean_disturbance_interval(30.0);
}

Disturbance::Disturbance(double mean_disturbance_interval)
  : shape(2.0) {
  set_mean_disturbance_interval(mean_disturbance_interval);
}

double Disturbance::survival_probability(double time_start,
					 double time) const {
  return survival0(time) / survival0(time_start);
}

double Disturbance::survival0(double time) const {
  return exp(-scale * pow(time, shape));
}

void Disturbance::set_mean_disturbance_interval(double x) {
  scale = pow(R::gammafn(1.0/shape)/shape/x, shape);
  p0 = shape*pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
}

}
