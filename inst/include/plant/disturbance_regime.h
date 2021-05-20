// -*-c++-*-
#ifndef PLANT_PLANT_DISTURBANCE_REGIME_H_
#define PLANT_PLANT_DISTURBANCE_REGIME_H_

#include <vector>
#include <RcppCommon.h> // NA_REAL

namespace plant {

class Disturbance_Regime {
public:
  virtual double density(double time) const {return NA_REAL;}
  virtual double pr_survival(double time) const {return NA_REAL;}
  virtual double cdf(double time) const {return NA_REAL;}
  virtual double icdf(double p) const {return NA_REAL;}
  virtual double r_mean_interval() const {return NA_REAL;}

  std::vector<double> r_density(std::vector<double> time) const;
};

}

#endif
