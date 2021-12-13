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

  virtual std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "net_mass_production_dt",
      "evapotranspiration"
    });
    // add the associated computation to compute_rates and compute there
    if (collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  double evapotranspiration_dt(double area_leaf) const;
  double water_access(const FF16_Environment& environment,
                      double height, double area_leaf_);

  // temporarily treating all soil layers as the same
  double compute_consumption(int resource_index) {
    return aux_index.at("evapotranspiration");
  }

  // Overrides ----------------------------------------------
  virtual double net_mass_production_dt(const FF16_Environment& environment,
                                        double height, double area_leaf_,
                                        bool reuse_intervals=false);



  virtual void compute_rates(const FF16_Environment& environment,
                             bool reuse_intervals,
                             Internals& vars);

};

FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s);

}

#endif
