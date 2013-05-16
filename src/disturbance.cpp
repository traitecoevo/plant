#include "disturbance.h"

#include <Rcpp.h>

namespace model {

// This should be rewritten to more closely use R's actual Weibull
// functions.

Disturbance::Disturbance(double mean_disturbance_interval) 
  : shape(2.0),
    scale(pow(R::gammafn(1.0/shape)/shape/mean_disturbance_interval, shape)),
    p0(shape*pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape)) {
}

double Disturbance::survival_probability(double time_start,
					 double time) const {
  return survival0(time) / survival0(time_start);
}

double Disturbance::survival0(double time) const {
  return exp(-scale * pow(time, shape));
}

}
