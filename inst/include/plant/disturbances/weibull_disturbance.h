// -*-c++-*-
#ifndef PLANT_PLANT_WEIBULL_DisturbanceRegime_H_
#define PLANT_PLANT_WEIBULL_DisturbanceRegime_H_

#include <vector>
#include <plant/disturbance_regime.h>
#include <Rcpp.h>

using namespace Rcpp;

namespace plant {

class Weibull_DisturbanceRegime: public DisturbanceRegime {
public:

  // icdf_limit chosen to match Falster 2011
  Weibull_DisturbanceRegime(double patch_max_lifetime)
    : DisturbanceRegime(),
      icdf_limit(6.25302620663814e-05),
      shape(2.0) {
    mean_interval = (patch_max_lifetime * R::gammafn(1.0 / shape)) /
                      (shape * sqrt(-log(icdf_limit)));
    scale = pow(R::gammafn(1 + 1.0 / shape) / mean_interval, shape);
    p0 = shape * pow(scale, 1.0 / shape) / R::gammafn(1.0 / shape);
  }

  // probability that a patch survives from time 0 to time
  virtual double pr_survival(double time) const {
    return exp(-scale * pow(time, shape));
  }

  // unnormalised patch density
  virtual double density(double time) const {
    return p0 * pr_survival(time);
  }

  // cumulative density function
  double cdf(double time) const {
    return 1 - pr_survival(time);
  }

  // Inverse cumulative density function
  double icdf(double p) const {
    return pow(log(p) / -scale, 1/shape);
  }

  // R accessor for distribution mean
  double r_mean_interval() const {
    return mean_interval;
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
