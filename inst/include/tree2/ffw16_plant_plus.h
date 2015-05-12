// -*-c++-*-
#ifndef TREE_FFW16_PLANT_PLUS_H_
#define TREE_FFW16_PLANT_PLUS_H_

#include <vector>
#include <tree2/ffw16_strategy.h>
#include <tree2/ffw16_plant_plus_internals.h>
#include <tree2/ode_interface.h>

namespace tree2 {

class FFW16_PlantPlus {
public:
  typedef FFW16_PlantPlus_internals internals;
  typedef FFW16_Strategy  strategy_type;
  FFW16_PlantPlus(strategy_type::ptr s);

  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.

  // * Individual size
  double height() const;
  double height_dt() const;
  void set_height(double height_);

  double mortality() const;
  double mortality_dt() const;
  void set_mortality(double x);

  double fecundity() const;
  double fecundity_dt() const;
  void set_fecundity(double x);

  double area_heartwood() const;
  double area_heartwood_dt() const;
  void set_area_heartwood(double x);

  double mass_heartwood() const;
  double mass_heartwood_dt() const;
  void set_mass_heartwood(double x);

  // * Competitive environment
  double area_leaf() const;
  // [      ] Leaf area (not fraction) above height `z`
  double area_leaf_above(double z) const;

  // * Mass production
  // [eqn 12-19,21] Update physiological variables
  void compute_vars_phys(const Environment& environment, bool
			 reuse_intervals=false);
  // TODO: This is temporary -- it should be called by
  // compute_vars_phys, but I don't want that to always happen, so do
  // it manually for now.
  void compute_vars_growth();

  // * Births and deaths
  // [eqn 20] Survival of seedlings during germination
  double germination_probability(const Environment& environment);

  // * ODE interface
  static size_t       ode_size() {return 5;}
  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;
  // Optional, but useful
  static std::vector<std::string> ode_names();

  // * R interface
  strategy_type r_get_strategy() const;

  // * Used by tools:
  double net_mass_production_dt() const {return vars.net_mass_production_dt;}

  FFW16_PlantPlus r_copy() const;
  internals r_internals() const;

  const Control& control() const;

private:
  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void compute_vars_size(double height_);

  strategy_type::ptr strategy;
  internals vars;

  // The ode dimensions are:
  // 1. Height
  // 2. Mortality
  // 3. Fecundity
  // 4. Heartwood area
  // 5. Heartwood mass
  static const int ode_dimension = 5;
};

FFW16_PlantPlus make_plant_plus(FFW16_PlantPlus::strategy_type s);

}

#endif
