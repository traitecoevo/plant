// Built from  inst/include/plant/models/ff16_strategy.h using the scaffolder,
// from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_WATER_STRATEGY_H_
#define PLANT_PLANT_WATER_STRATEGY_H_

#include <plant/leaf_model.h>
#include <plant/models/ff16_strategy.h>

namespace plant {

class FF16w_Strategy : public FF16_Strategy {
public:
  typedef std::shared_ptr<FF16w_Strategy> ptr;
  FF16w_Strategy();

  double compute_average_light_environment(double z, double height,
                                           const FF16_Environment &environment);

  double evapotranspiration_dt(double area_leaf_);

  virtual double net_mass_production_dt(const FF16_Environment &environment,
                                        double height, double area_leaf_,
                                        bool reuse_intervals = false);


  virtual void compute_rates(const FF16_Environment &environment,
                             bool reuse_intervals, Internals &vars);

  Leaf leaf;

  // leaf traits
  double vcmax = 100;
  double p_50 = 1.0;
  double c = 2.04;
  double b = 2.0;
  double psi_crit = 3.42; // derived from b and c
  double beta = 15000.0;
  double beta_2 = 1.0;
  double huber_value = 1.57e-4;
};

FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s);

} // namespace plant

#endif
