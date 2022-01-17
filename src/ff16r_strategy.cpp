// Built from  src/ff16_strategy.cpp on Fri Jul  3 08:14:35 2020 using the scaffolder, from the strategy:  FF16
#include <plant/models/ff16r_strategy.h>

namespace plant {

FF16r_Strategy::FF16r_Strategy() {
  collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16r";
}

// probabilty of establishment decays over time
double FF16r_Strategy::establishment_probability(const FF16_Environment& environment) {
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0);

  double decay_over_time = exp(-recruitment_decay * environment.time);

  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return (1.0 / (tmp * tmp + 1.0)) * decay_over_time;
  } else {
    return 0.0;
  }
}


FF16r_Strategy::ptr make_strategy_ptr(FF16r_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16r_Strategy>(s);
}
}
