// -*-c++-*-
#ifndef PLANT_PLANT_NO_DISTURBANCE_H_
#define PLANT_PLANT_NO_DISTURBANCE_H_

#include <vector>
#include <plant/disturbance_regime.h>

namespace plant {

class No_Disturbance: public Disturbance_Regime {
public:
  No_Disturbance();
  virtual double density(double time) const;
  virtual double pr_survival(double time) const;
};

}

#endif
