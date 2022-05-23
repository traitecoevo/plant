// Built from  inst/include/plant/models/ff16_strategy.h on Fri Oct 30 11:43:30
// 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_FF16BG_STRATEGY_H_
#define PLANT_PLANT_FF16BG_STRATEGY_H_

#include <plant/models/assimilation.h>
#include <plant/models/ff16_environment.h>
#include <plant/models/ff16_strategy.h>

namespace plant {

class FF16bg_Strategy : public FF16_Strategy {
public:
  typedef std::shared_ptr<FF16bg_Strategy> ptr;
  FF16bg_Strategy();

  // Overrides ----------------------------------------------
  // Controls the relative effect of belowground competition
  double k_2 = 0;

  // extrinsic mortality to complement intrinsic and growth related mortality
  // default is constant 0.0
  std::vector<double> mortality_rate_x;
  std::vector<double> mortality_rate_y = {0.0};
  // whether the spline for each species should be constant fn or not
  // (extrapolation on/off)
  bool is_variable_mortality_rate = false;

  virtual double mortality_dt(double productivity_area,
                              double cumulative_mortality, double time) const;

  double below_ground_influence(const FF16_Environment &environment);

  // Incorporates FF16 extensions above
  virtual double net_mass_production_dt(const FF16_Environment &environment,
                                        double height, double area_leaf_,
                                        bool reuse_intervals = false);

  // one-shot update of the scm variables
  // i.e. setting rates of ode vars from the state and updating aux vars
  virtual void compute_rates(const FF16_Environment &environment,
                             bool reuse_intervals, Internals &vars);

  void prepare_strategy();
};

FF16bg_Strategy::ptr make_strategy_ptr(FF16bg_Strategy s);

} // namespace plant

#endif
