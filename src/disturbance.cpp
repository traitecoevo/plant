#include "disturbance.h"

#include <Rcpp.h>

namespace model {

Disturbance::Disturbance(double mean_interval_)
  : shape(2.0),
    mean_interval(mean_interval_) {
  scale = pow(R::gammafn(1.0/shape)/shape/mean_interval, shape);
  p0 = shape*pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
}

double Disturbance::survival_probability(double time_start,
					 double time) const {
  return survival0(time) / survival0(time_start);
}

double Disturbance::r_mean_interval() const {
  return mean_interval;
}

double Disturbance::survival0(double time) const {
  return exp(-scale * pow(time, shape));
}

}
