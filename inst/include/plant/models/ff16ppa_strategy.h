// Built from  inst/include/plant/models/ff16_strategy.h on Thu Oct 29 11:14:41 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_FF16ppa_STRATEGY_H_
#define PLANT_PLANT_FF16ppa_STRATEGY_H_

#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>

namespace plant {

class FF16ppa_Strategy: public FF16_Strategy {
public:
  typedef std::shared_ptr<FF16ppa_Strategy> ptr;
  FF16ppa_Strategy();

  // Overrides ----------------------------------------------
  // Perfect plasticity approximation
  virtual double ppa(const FF16_Environment& environment,
                      double height, double area_leaf_);

  // [eqn 16] Fraction of whole plant growth that is leaf
  virtual double net_mass_production_dt(const FF16_Environment& environment,
                                        double height, double area_leaf_,
                                        bool reuse_intervals=false);
};

FF16ppa_Strategy::ptr make_strategy_ptr(FF16ppa_Strategy s);

}

#endif
