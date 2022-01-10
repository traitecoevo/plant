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

  double evapotranspiration_dt(double area_leaf) const;

  virtual void compute_rates(const FF16_Environment& environment,
                             bool reuse_intervals,
                             Internals& vars);
};

FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s);

}

#endif
