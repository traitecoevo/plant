// Built from  src/ff16_strategy.cpp on Fri Oct 30 11:43:30 2020 using the scaffolder, from the strategy:  FF16
#include <plant/models/ff16bg_strategy.h>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
FF16bg_Strategy::FF16bg_Strategy() {
  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16bg";
}

double FF16bg_Strategy::below_ground_influence(const FF16_Environment& environment) {
  double e_0 = environment.get_environment_at_height(height_0 + 1e-6);
  return(pow(e_0, k_2 / environment.canopy.k_I));
}

double FF16bg_Strategy::establishment_probability(const FF16_Environment& environment) {
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0);

  double x = exp(-a_f4 * environment.time);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * x;
  } else {
    return 0.0;
  }
}

// Add rate limit based on total stand leaf area
double FF16bg_Strategy::net_mass_production_dt(const FF16_Environment& environment,
                                double height, double area_leaf_,
                                bool reuse_intervals) {
  const double mass_leaf_    = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);

  // Reduce net assimilation of all individuals
  const double assimilation_ = assimilator.assimilate(control, environment, height,
                                            area_leaf_, reuse_intervals) *
                                below_ground_influence(environment);
  const double respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
    turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

FF16bg_Strategy::ptr make_strategy_ptr(FF16bg_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16bg_Strategy>(s);
}
}
