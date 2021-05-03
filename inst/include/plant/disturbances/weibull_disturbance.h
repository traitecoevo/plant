// -*-c++-*-
#ifndef PLANT_PLANT_WEIBULL_DISTURBANCE_REGIME_H_
#define PLANT_PLANT_WEIBULL_DISTURBANCE_REGIME_H_

#include <vector>
#include <plant/disturbance_regime.h>
#include <Rcpp.h>

using namespace Rcpp;

namespace plant {

class Weibull_Disturbance_Regime: public Disturbance_Regime {
public:
  Weibull_Disturbance_Regime()
    : shape(2.0),
      mean_interval(1.0) {
    scale = pow(R::gammafn(1.0/shape)/shape/mean_interval, shape);
    p0 = shape*pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
  }

  Weibull_Disturbance_Regime(double mean_interval_)
    : shape(2.0),
      mean_interval(mean_interval_) {
    scale = pow(R::gammafn(1.0/shape)/shape/mean_interval, shape);
    p0 = shape*pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
  }

  virtual double density(double time) const {
    return p0 * pr_survival(time);
  }

  // probability that a patch survives from time 0 to time
  virtual double pr_survival(double time) const {
    return exp(-scale * pow(time, shape));
  }

  double r_mean_interval() const {
    return mean_interval;
  }

  // Inverse cumulative density function
  double icdf(double p) const {
    return pow(log(p) / -scale, 1/shape);
  }

private:
  double shape;
  double mean_interval;
  double scale;
  double p0;
};

}

#endif
