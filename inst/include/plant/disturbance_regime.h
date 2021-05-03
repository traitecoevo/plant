// -*-c++-*-
#ifndef PLANT_PLANT_DISTURBANCE_REGIME_H_
#define PLANT_PLANT_DISTURBANCE_REGIME_H_

#include <vector>

namespace plant {

class Disturbance_Regime {
public:
  double density(double time) const;
  double pr_survival(double time) const;
  std::vector<double> r_density(std::vector<double> time) const;
};

}

#endif
