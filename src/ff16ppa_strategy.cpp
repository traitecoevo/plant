// Built from  src/ff16_strategy.cpp on Thu Oct 29 11:14:41 2020 using the scaffolder, from the strategy:  FF16
#include <plant/models/ff16ppa_strategy.h>

namespace plant {

FF16ppa_Strategy::FF16ppa_Strategy() {
  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16ppa";
}

// Introduce step function
double FF16ppa_Strategy::ppa(const FF16_Environment& environment,
                    double height, double area_leaf) {

  // spline sometimes > 1
  double E = std::min(environment.get_environment_at_height(height), 1.0);

  // back-transform and truncate
  double L = - log(E) / 0.5;
  double E2 = exp(- 0.5 * floor(L));

  return (a_p1 * E2 / (E2 + a_p2)) * area_leaf;
}

// Use PPA for net assimilation
double FF16ppa_Strategy::net_mass_production_dt(const FF16_Environment& environment,
                                  double height, double area_leaf_, bool reuse_intervals) {
  const double mass_leaf_    = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);

  // replace assimilator
  const double assimilation_ = ppa(environment, height, area_leaf_);

  const double respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
    turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

FF16ppa_Strategy::ptr make_strategy_ptr(FF16ppa_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16ppa_Strategy>(s);
}
}
