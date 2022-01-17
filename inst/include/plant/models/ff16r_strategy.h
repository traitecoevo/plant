// Built from  inst/include/plant/models/ff16_strategy.h on Fri Jul  3 08:14:35 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_FF16R_STRATEGY_H_
#define PLANT_PLANT_FF16R_STRATEGY_H_

#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/models/assimilation.h>

namespace plant {

class FF16r_Strategy: public FF16_Strategy {
public:
  typedef std::shared_ptr<FF16r_Strategy> ptr;

  FF16r_Strategy();

  // Overloads ----------------------------------------------

  // [eqn 20] Survival of seedlings during establishment
  virtual double establishment_probability(const FF16_Environment& environment);

  double recruitment_decay = 1;

};

FF16r_Strategy::ptr make_strategy_ptr(FF16r_Strategy s);

}

#endif
