// Built from  inst/include/plant/models/ff16_strategy.h using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_WATER_STRATEGY_H_
#define PLANT_PLANT_WATER_STRATEGY_H_

#include <plant/models/ff16_strategy.h>

namespace plant {

class FF16w_Strategy: public FF16_Strategy {
public:
  typedef std::shared_ptr<FF16w_Strategy> ptr;
  FF16w_Strategy();

  // Overrides ----------------------------------------------
  virtual double net_mass_production_dt(const FF16_Environment& environment,
                                        double height, double area_leaf_,
                                        bool reuse_intervals=false);

  double water_access(const FF16_Environment& environment,
                      double height, double area_leaf_);
};

FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s);

}

#endif
