// -*-c++-*-
#ifndef TREE_PLANT_H_
#define TREE_PLANT_H_

#include <memory> // std::shared_ptr
#include <vector>
#include <tree2/strategy.h>
#include <tree2/plant_internals.h>
#include <tree2/ode_interface.h>
#include <RcppCommon.h>

namespace tree2 {

class Plant {
public:
  typedef Strategy                       strategy_type;
  typedef std::shared_ptr<strategy_type> strategy_ptr_type;
  Plant(strategy_ptr_type s);

  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.

  // * Individual size
  double height() const;
  double height_rate() const;
  void set_height(double height_);

  double mortality() const;
  double mortality_rate() const;
  void set_mortality(double x);

  double fecundity() const;
  double fecundity_rate() const;
  void set_fecundity(double x);

  double heartwood_area() const;
  double heartwood_area_rate() const;
  void set_heartwood_area(double x);

  double heartwood_mass() const;
  double heartwood_mass_rate() const;
  void set_heartwood_mass(double x);

  // These are derived from mortality() -- see design.md.
  double survival_probability() const;

  // * Competitive environment
  double leaf_area() const;
  // [      ] Leaf area (not fraction) above height `z`
  double leaf_area_above(double z) const;

  // * Mass production
  // [eqn 12-19,21] Update physiological variables
  void compute_vars_phys(const Environment& environment, bool
			 reuse_intervals=false);

  // * Births and deaths
  // [eqn 20] Survival of seedlings during germination
  double germination_probability(const Environment& environment);

  // * ODE interface
  static size_t ode_size() {return 5;}
  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;
  // Optional, but useful
  static std::vector<std::string> ode_names();

  // * R interface
  strategy_type r_get_strategy() const;
  SEXP r_get_vars_growth() const;

  // * Used by tools:
  double net_production() const {return vars.net_production;}

  Plant r_copy() const;
  Plant_internals r_internals() const;

  const Control& control() const;

private:
  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void compute_vars_size(double height_);

  strategy_ptr_type strategy;
  Plant_internals vars;

  // The ode dimensions are:
  // 1. Height
  // 2. Mortality
  // 3. Fecundity
  // 4. Heartwood area
  // 5. Heartwood mass
  static const int ode_dimension = 5;
};

Plant::strategy_ptr_type make_strategy_ptr(Plant::strategy_type s);
Plant make_plant(Plant::strategy_type s);

}

#endif
