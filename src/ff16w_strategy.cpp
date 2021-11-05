// Inherit from FF16, like FF16r
#include <plant/models/ff16w_strategy.h>

namespace plant {
FF16w_Strategy::FF16w_Strategy() {
  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16w";
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
double FF16w_Strategy::net_mass_production_dt(const FF16_Environment& environment,
                                double height, double area_leaf_,
                                bool reuse_intervals) {
  const double mass_leaf_    = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);
  const double assimilation_ = assimilator.assimilate(environment, height,
                                                      area_leaf_, reuse_intervals);
  const double respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
    turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);

  const double water_ = water_access(environment, height, area_leaf_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// TODO: basic test of accessing environment state
double FF16w_Strategy::water_access(const FF16_Environment &environment,
                                    double height, double area_leaf_) {

  std::vector<double> water_ = environment.get_soil_water_state();
  return water_[0];
}

FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16w_Strategy>(s);
}
}
