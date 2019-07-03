// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_MINIMAL_H_
#define PLANT_PLANT_PLANT_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <plant/ode_interface.h>
#include <vector>
#include <plant/internals.h>

namespace plant {

template <typename T> class Plant {
public:
  typedef T strategy_type;
  typedef Internals internals;
  typedef typename strategy_type::ptr strategy_type_ptr;
  // for the time being...
  Plant(strategy_type_ptr s) : strategy(s) {
    set_state("height", strategy->height_0);
  }

  // useage: state("height")
  double state(std::string name) const {
    return vars.state(strategy->state_index[name]);
  }
  
  // useage:_rate("area_heartwood")
  double rate(std::string name) const {
    return vars.rate(strategy->state_index[name]);
  }

  // useage: set_state("height", 2.0)
  void set_state(std::string name, double v) {
    vars.set_state(strategy->state_index[name], v);
  }

  double area_leaf_above(double z) const {
    return strategy->area_leaf_above(z, state("height"));
  }

  void compute_vars_phys(const Environment &environment,
                         bool reuse_intervals = false) {
    strategy->scm_vars(environment, reuse_intervals, vars);
  }
  double germination_probability(const Environment &environment) {
    return strategy->germination_probability(environment);
  }

  // * ODE interface
  static size_t ode_size() { return 5; } // we want it to be: strategy->state_size; }

  ode::const_iterator set_ode_state(ode::const_iterator it) {
    for (int i = 0; i < vars.state_size; i++) {
      vars.states[i] = *it++;
    }
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    for (int i = 0; i < vars.state_size; i++) {
      *it++ = vars.states[i];
    }
    return it;
  }
  ode::iterator ode_rates(ode::iterator it) const {
    for (int i = 0; i < vars.state_size; i++) {
      *it++ = vars.rates[i];
    }
    return it;
  }
  // Optional, but useful
  static std::vector<std::string> ode_names() { return strategy_type::state_names(); }

  // Used in the stochastic model:
  double mortality_probability() const { return 1 - exp(-state("mortality")); }
  
  void reset_mortality() { set_state("mortality", 0.0); }

  std::string strategy_name() const { return strategy->name; }

  // * R interface
  strategy_type r_get_strategy() const { return *strategy.get(); }
  // ! External R code depends on knowing r internals for like growing plant to
  // ! height or something
  internals r_internals() const { return vars; }
  const Control &control() const { return strategy->control; }

private:
  strategy_type_ptr strategy;
  internals vars;
};

template <typename T> Plant<T> make_plant(T s) {
  return Plant<T>(make_strategy_ptr(s));
}

} // namespace plant

#endif
