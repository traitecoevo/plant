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

  // icdf_limit chosen to match Falster 2011
  Weibull_Disturbance_Regime(double max_patch_lifetime)
    : icdf_limit(6.25302620663814e-05),
      shape(2.0) {
    mean_interval = (max_patch_lifetime * R::gammafn(1.0 / shape)) /
                      (shape * sqrt(-log(icdf_limit)));
    scale = pow(R::gammafn(1.0 / shape) / shape / mean_interval, shape);
    p0 = shape * pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
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
  double icdf_limit;
  double mean_interval;
  double shape;
  double scale;
  double p0;
};

}

#endif
