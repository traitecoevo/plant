#include "disturbance.h"

#include <Rcpp.h>

namespace model {

Disturbance::Disturbance()
  : shape(2.0),
    mean_disturbance_interval(30) {
  set_parameters_post_hook();
}

double Disturbance::survival_probability(double time_start,
					 double time) const {
  return survival0(time) / survival0(time_start);
}

double Disturbance::survival0(double time) const {
  return exp(-scale * pow(time, shape));
}

void Disturbance::do_build_lookup() {
  lookup_table["mean_disturbance_interval"] = &mean_disturbance_interval;
}

void Disturbance::set_parameters_post_hook() {
  scale = pow(R::gammafn(1.0/shape)/shape/mean_disturbance_interval, shape);
  p0 = shape*pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
}

}
