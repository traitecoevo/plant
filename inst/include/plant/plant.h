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
  typedef typename strategy_type::ptr strategy_type_ptr;
  // for the time being...
  Plant(strategy_type_ptr s) : strategy(s) {
    vars.resize(strategy_type::state_size(), strategy_type::aux_size()); // = Internals(strategy_type::state_size());
    set_state("height", strategy->height_0);
  }

  // useage: state(HEIGHT_INDEX)
  double state(std::string name) const {
    return vars.state(strategy->state_index.at(name));
  }
  double state(int i) const { return vars.state(i); }
  
  // useage:_rate("area_heartwood")
  double rate(std::string name) const {
    return vars.rate(strategy->state_index.at(name));
  }
  double rate(int i) const { return vars.rate(i); }

  // useage: set_state("height", 2.0)
  void set_state(std::string name, double v) {
    vars.set_state(strategy->state_index.at(name), v);
  }
  void set_state(int i, double v) { vars.set_state(i, v); }

  // aux vars by name and index
  double aux(std::string name) const {
    return vars.aux(strategy->aux_index.at(name));
  }
  double aux(int i) const { return vars.aux(i); } 

  double area_leaf_above(double z) const {
    return strategy->area_leaf_above(z, state(HEIGHT_INDEX));
  }

  void compute_vars_phys(const Environment &environment,
                         bool reuse_intervals = false) {
    strategy->compute_vars_phys(environment, reuse_intervals, vars);
  }
  double germination_probability(const Environment &environment) {
    return strategy->germination_probability(environment);
  }

  // * ODE interface
  static size_t ode_size() { return strategy_type::state_size(); }
  static size_t aux_size() { return strategy_type::aux_size(); }
  static std::vector<std::string> aux_names() { return strategy_type::aux_names(); }
  static std::vector<std::string> ode_names() { return strategy_type::state_names(); }

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

  // Used in the stochastic model:
  double mortality_probability() const { return 1 - exp(-state(MORTALITY_INDEX)); }
  
  void reset_mortality() { set_state("mortality", 0.0); }

  std::string strategy_name() const { return strategy->name; }

  // * R interface
  strategy_type r_get_strategy() const { return *strategy.get(); }
  // ! External R code depends on knowing r internals for like growing plant to
  // ! height or something
  Internals r_internals() const { return vars; }
  const Control &control() const { return strategy->control; }

private:
  strategy_type_ptr strategy;
  Internals vars;
};

template <typename T> Plant<T> make_plant(T s) {
  return Plant<T>(make_strategy_ptr(s));
}

} // namespace plant

#endif
