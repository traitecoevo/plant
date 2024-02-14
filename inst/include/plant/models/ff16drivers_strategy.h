// -*-c++-*-
#ifndef PLANT_PLANT_FF16DRIVERS_STRATEGY_H_
#define PLANT_PLANT_FF16DRIVERS_STRATEGY_H_

#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/models/assimilation.h>

namespace plant {

class FF16drivers_Strategy: public FF16_Strategy {
public:
  typedef std::shared_ptr<FF16drivers_Strategy> ptr;
  FF16drivers_Strategy();

  // add external drivers
  std::vector<double> growth_rate_x;
  std::vector<double> growth_rate_y = {1.0};
  bool is_variable_growth_rate = false;

  std::vector<double> mortality_rate_x;
  std::vector<double> mortality_rate_y = {1.0};
  bool is_variable_mortality_rate = false;

  // Overloads ----------------------------------------------
  virtual void compute_rates(const FF16_Environment& environment, bool reuse_intervals,
                              Internals& vars);

  virtual double net_mass_production_dt(const FF16_Environment& environment,
                                        double height, double area_leaf_, double time,
                                        bool reuse_intervals=false);

  virtual double mortality_dt(double productivity_area, double cumulative_mortality, double time, double height) const;
  virtual double mortality_growth_independent_dt(double time, double height) const;
  
  virtual double establishment_probability(const FF16_Environment& environment);

  void prepare_strategy();
};

FF16drivers_Strategy::ptr make_strategy_ptr(FF16drivers_Strategy s);

}

#endif
