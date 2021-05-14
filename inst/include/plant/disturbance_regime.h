// -*-c++-*-
#ifndef PLANT_PLANT_DISTURBANCE_REGIME_H_
#define PLANT_PLANT_DISTURBANCE_REGIME_H_

#include <vector>
#include <RcppCommon.h> // NA_REAL

namespace plant {

class Disturbance_Regime {
public:
  virtual double density(double time) const {return NA_REAL;};
  virtual double pr_survival(double time) const {return NA_REAL;};
  std::vector<double> r_density(std::vector<double> time) const;

};

}

#endif
