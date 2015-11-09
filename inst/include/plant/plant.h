// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_MINIMAL_H_
#define PLANT_PLANT_PLANT_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <vector>
#include <plant/ode_interface.h>
#include <plant/plant_internals.h>

namespace plant {

template <typename T>
class Plant {
public:
  typedef Plant_internals internals;
  typedef T               strategy_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Plant(strategy_type_ptr s)
    : strategy(s) {
    set_height(strategy->height_0);
  }

  double height() const {return vars.height;}
  double height_dt() const {return vars.height_dt;}
  void set_height(double x) {
    vars.height    = x;
    vars.area_leaf = strategy->area_leaf(x);
  }

  double mortality() const {return vars.mortality;}
  double mortality_dt() const {return vars.mortality_dt;}
  void set_mortality(double x) {vars.mortality = x;}

  double fecundity() const {return vars.fecundity;}
  double fecundity_dt() const {return vars.fecundity_dt;}
  void set_fecundity(double x) {vars.fecundity = x;}

  double area_heartwood() const {return vars.area_heartwood;}
  double area_heartwood_dt() const {return vars.area_heartwood_dt;}
  void set_area_heartwood(double x) {vars.area_heartwood = x;}

  double mass_heartwood() const {return vars.mass_heartwood;}
  double mass_heartwood_dt() const {return vars.mass_heartwood_dt;}
  void set_mass_heartwood(double x) {vars.mass_heartwood = x;}

  double area_leaf_above(double z) const {
    return strategy->area_leaf_above(z, vars.height, vars.area_leaf);
  }

  void compute_vars_phys(const Environment& environment,
                         bool reuse_intervals=false) {
    strategy->scm_vars(environment, reuse_intervals, vars);
  }
  double germination_probability(const Environment& environment) {
    return strategy->germination_probability(environment);
  }

  // * ODE interface
  static size_t       ode_size() {return 5;}
  ode::const_iterator set_ode_state(ode::const_iterator it) {
    set_height(*it++);
    set_mortality(*it++);
    set_fecundity(*it++);
    set_area_heartwood(*it++);
    set_mass_heartwood(*it++);
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    *it++ = height();
    *it++ = mortality();
    *it++ = fecundity();
    *it++ = area_heartwood();
    *it++ = mass_heartwood();
    return it;
  }
  ode::iterator ode_rates(ode::iterator it) const {
    *it++ = height_dt();
    *it++ = mortality_dt();
    *it++ = fecundity_dt();
    *it++ = area_heartwood_dt();
    *it++ = mass_heartwood_dt();
    return it;
  }
  // Optional, but useful
  static std::vector<std::string> ode_names() {
    return std::vector<std::string>({"height", "mortality", "fecundity",
          "area_heartwood", "mass_heartwood"});
  }

  // Used in the stochastic model:
  double mortality_probability() const {
    return 1 - exp(-mortality());
  }
  void reset_mortality() {
    set_mortality(0.0);
  }

  // * R interface
  strategy_type r_get_strategy() const {return *strategy.get();}
  internals r_internals() const {return vars;}
  const Control& control() const {return strategy->control;}

private:
  strategy_type_ptr strategy;
  internals vars;
};

template <typename T>
Plant<T> make_plant(T s) {
  return Plant<T>(make_strategy_ptr(s));
}

}

#endif
