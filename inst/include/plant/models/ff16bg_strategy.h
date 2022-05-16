// Built from  inst/include/plant/models/ff16_strategy.h on Fri Oct 30 11:43:30 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_FF16BG_STRATEGY_H_
#define PLANT_PLANT_FF16BG_STRATEGY_H_

#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/models/assimilation.h>

namespace plant {

class FF16bg_Strategy: public FF16_Strategy {
public:
  typedef std::shared_ptr<FF16bg_Strategy> ptr;
  FF16bg_Strategy();

  // Overrides ----------------------------------------------
  // Controls the relative effect of belowground competition
  double k_2 = 0;

  double below_ground_influence(const FF16_Environment& environment);

  // Incorporates FF16 extensions above
  virtual double net_mass_production_dt(const FF16_Environment& environment,
                                        double height, double area_leaf_,
                                        bool reuse_intervals=false);
};

FF16bg_Strategy::ptr make_strategy_ptr(FF16bg_Strategy s);

}

#endif
