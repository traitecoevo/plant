// -*-c++-*-
#ifndef PLANT_PLANT_STRATEGY_H_
#define PLANT_PLANT_STRATEGY_H_

#include <memory>
#include <plant/control.h>
#include <plant/qag_internals.h> // quadrature::intervals_type
#include <plant/internals.h> // quadrature::intervals_type
#include <plant/strategy.h>

namespace plant {

class Kohyama_Strategy: public Strategy {
public:
  typedef std::shared_ptr<Strategy> ptr;
  Kohyama_Strategy();

  double competition_effect(double size) const;

  double competition_effect_state(Internals& vars);

  void compute_rates(const Environment& environment, bool reuse_intervals,
                Internals& vars);

  double net_mass_production_dt(const Environment& environment,
                                double size, double competition_effect_,
                                bool reuse_intervals=false);

  double germination_probability(const Environment& environment);

  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  double mortality_dt(double productivity_area, double cumulative_mortality) const;

  double compute_competition(double z, double size, double competition_effect) const;

  double basal_area(double height) {

  }

};

Strategy::ptr make_strategy_ptr(Strategy s);

}

#endif
