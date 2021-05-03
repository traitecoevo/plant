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
  Weibull_Disturbance_Regime();
  Weibull_Disturbance_Regime(double mean_interval_);

  // overloaded from base class
  virtual double density(double time) const;
  virtual double pr_survival(double time) const;

  double r_mean_interval() const;
  double icdf(double pr) const;

private:
  double shape;
  double mean_interval;
  double scale;
  double p0;
};

}

#endif
