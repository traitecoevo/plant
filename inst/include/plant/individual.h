// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_MINIMAL_H_
#define PLANT_PLANT_PLANT_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <odelia/ode_interface.hpp>
#include <vector>
#include <plant/internals.h>
#include <plant/uniroot.h>


namespace plant {

template <typename T, typename E> class Individual {
public:
  typedef T strategy_type;
  typedef E environment_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  // for the time being...
  Individual(strategy_type_ptr s) : strategy(s) {
    if (strategy->aux_index.size() != s->aux_size()) {
      strategy->refresh_indices();
    }
    // Resolve the named aux slots once at construction so the hot
    // compute_competition() / net_mass_production_dt() paths read them by
    // integer index instead of a std::map<string,int>::at lookup per call
    // (those lookups were visible in profiling, see #466).
    competition_effect_aux_index = strategy->aux_index.at("competition_effect");
    height_inverse_aux_index = strategy->aux_index.at("height_inverse");
    vars.resize(strategy_type::state_size(), s->aux_size()); // = Internals(strategy_type::state_size());
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
    int i = strategy->state_index.at(name);
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }
  void set_state(int i, double v) {
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }

  // aux vars by name and index
  double aux(std::string name) const {
    return vars.aux(strategy->aux_index.at(name));
  }
  double aux(int i) const { return vars.aux(i); }

  // set # consumable resources based on env. variables
  void resize_consumption_rates(int i) {
    vars.resize_consumption_rates(i);
  }
  double consumption_rate(int i) const { return vars.consumption_rate(i); }

  double compute_competition(double z) const {
    return strategy->compute_competition(
      z,
      vars.aux(competition_effect_aux_index),
      vars.aux(height_inverse_aux_index));
  }

  void compute_rates(const environment_type& environment) {
    if (vars.resource_size != environment.ode_size()) {
      // handles when Individual hasn't been instantiated in a Patch (ie with an environment)
      vars.resize_consumption_rates(environment.ode_size());
    }
    strategy->compute_rates(environment, vars);
  }
  
  double establishment_probability(const environment_type &environment) {
    return strategy->establishment_probability(environment);
  }

  double net_mass_production_dt(const environment_type &environment) {
    // TODO(#483):  maybe reuse intervals? default false 
    return strategy->net_mass_production_dt(
      environment,
      state(HEIGHT_INDEX),
      vars.aux(competition_effect_aux_index),
      vars.aux(height_inverse_aux_index));
  }

  // * ODE interface
  static size_t ode_size() { return strategy_type::state_size(); }
  static std::vector<std::string> ode_names() { return strategy_type::state_names(); }

  // why doesn't strategy->aux_size() also need to be const-qualified? is it because strategy is a pointer?
  size_t aux_size() const { return strategy->aux_size(); }
  std::vector<std::string> aux_names() { return strategy->aux_names(); }

  odelia::ode::const_iterator set_ode_state(odelia::ode::const_iterator it) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.states[i] = *it++;
      strategy->update_dependent_aux(i, vars);
    }
    return it;
  }
  odelia::ode::iterator ode_state(odelia::ode::iterator it) const {
    for (size_t i = 0; i < vars.state_size; i++) {
      *it++ = vars.states[i];
    }
    return it;
  }
  odelia::ode::iterator ode_rates(odelia::ode::iterator it) const {
    for (size_t i = 0; i < vars.state_size; i++) {
      *it++ = vars.rates[i];
    }
    return it;
  }

  odelia::ode::iterator ode_aux(odelia::ode::iterator it) const {
    for (size_t i = 0; i < vars.aux_size; i++) {
      *it++ = vars.auxs[i];
    }
    return it;
  }

  // Single individual methods

  // Used in the stochastic model:
  double mortality_probability() const { return 1 - exp(-state(MORTALITY_INDEX)); }
  
  void reset_mortality() { set_state("mortality", 0.0); }

  double growth_rate_given_height(double height, const environment_type& environment) {
    // Called repeatedly from the finite-difference gradient (Node::
    // growth_rate_gradient), so address height by integer slot rather than the
    // "height" string-map lookup (see #466).
    set_state(HEIGHT_INDEX, height);
    compute_rates(environment);
    return rate(HEIGHT_INDEX);
  }

  double resource_compensation_point() {
    environment_type env = environment_type();

    auto target = [&] (double x) mutable -> double {
      env.set_fixed_environment(x, 100);
      compute_rates(env);
      return net_mass_production_dt(env);
    };

    const double f1 = target(1.0);
    if (f1 < 0.0) {
      return NA_REAL;
    } else {
      const double tol = control().offspring_production_tol;
      const size_t max_iterations = control().offspring_production_iterations;
      return util::uniroot(target, 0.0, 1.0, tol, max_iterations);
    }
  }



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
  // Cached aux slot indices (see constructor) for hot-path access.
  int competition_effect_aux_index;
  int height_inverse_aux_index;
};

template <typename T, typename E> Individual<T,E> make_individual(T s) {
  return Individual<T,E>(make_strategy_ptr(s));
}

} // namespace plant

#endif
