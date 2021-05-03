#include <plant/disturbances/weibull_disturbance.h>

namespace plant {

Weibull_Disturbance_Regime::Weibull_Disturbance_Regime()
  : shape(2.0),
    mean_interval(1.0) {
  scale = pow(R::gammafn(1.0/shape)/shape/mean_interval, shape);
  p0 = shape*pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
}

Weibull_Disturbance_Regime::Weibull_Disturbance_Regime(double mean_interval_)
  : shape(2.0),
    mean_interval(mean_interval_) {
  scale = pow(R::gammafn(1.0/shape)/shape/mean_interval, shape);
  p0 = shape*pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
}

double Weibull_Disturbance_Regime::density(double time) const {
  return p0 * pr_survival(time);
}

double Weibull_Disturbance_Regime::r_mean_interval() const {
  return mean_interval;
}

// This is the probability that a patch survives from time 0 to time
// 'time'.
double Weibull_Disturbance_Regime::pr_survival(double time) const {
  return exp(-scale * pow(time, shape));
}

// Inverse cumulative density function
double Weibull_Disturbance_Regime::icdf(double p) const {
  return pow(log(p) / -scale, 1/shape);
}

}
