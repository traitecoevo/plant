// -*-c++-*-
#ifndef PLANT_PLANT_NO_DISTURBANCE_H_
#define PLANT_PLANT_NO_DISTURBANCE_H_

#include <vector>
#include <plant/disturbance_regime.h>

namespace plant {


class No_Disturbance: public DisturbanceRegime {
public:
  virtual double density(double time) const {return 1.0;};
  virtual double pr_survival(double time) const {return 1.0;};
};

}

#endif
