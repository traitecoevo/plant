// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_MINIMAL_H_
#define PLANT_PLANT_PLANT_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <odelia/ode_interface.hpp>
#include <vector>
#include <plant/internals.h>
#include <plant/uniroot.h>


namespace plant {

// Templated on the scalar type S (#472 scope B / #537, Milestone C) so a plant's
// ODE state can be an AD active type for reverse-mode gradients. S defaults to
// double, so every existing `Individual<T,E>` is unchanged and bit-identical;
// only AD paths instantiate Individual<T,E,ad_type> (and then only the members
// they use are compiled -- the double-only ODE-iterator methods stay uncompiled
// for the AD instantiation until the ODE-state boundary is wired).
template <typename T, typename E, typename S = double> class Individual {
public:
  typedef T strategy_type;
  typedef E environment_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  // for the time being...
  Individual(strategy_type_ptr s) : strategy(s) {
    if (strategy->aux_index.size() != s->aux_size()) {
      strategy->refresh_indices();
    }
    // The strategy owns the knowledge of which state/aux slots its rate
    // functions need: Individual passes the whole Internals to the strategy
    // (see compute_competition / net_mass_production_dt below) rather than
    // resolving strategy-specific aux slots here (#266). The strategy still
    // reads them by cached integer index, so the #466 hot-path optimisation is
    // preserved.
    vars.resize(strategy_type::state_size(), s->aux_size()); // = Internals(strategy_type::state_size());
    set_state("height", strategy->initial_height());
  }
  
  // useage: state(HEIGHT_INDEX)
  S state(std::string name) const {
    return vars.state(strategy->state_index.at(name));
  }
  S state(int i) const { return vars.state(i); }

  // useage:_rate("area_heartwood")
  S rate(std::string name) const {
    return vars.rate(strategy->state_index.at(name));
  }
  S rate(int i) const { return vars.rate(i); }

  // useage: set_state("height", 2.0)
  void set_state(std::string name, S v) {
    int i = strategy->state_index.at(name);
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }
  void set_state(int i, S v) {
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }

  // aux vars by name and index
  S aux(std::string name) const {
    return vars.aux(strategy->aux_index.at(name));
  }
  S aux(int i) const { return vars.aux(i); }

  // set # consumable resources based on env. variables
  void resize_consumption_rates(int i) {
    vars.resize_consumption_rates(i);
  }
  S consumption_rate(int i) const { return vars.consumption_rate(i); }

  double compute_competition(double z) const {
    return strategy->compute_competition(z, vars);
  }

  // Seed strategy-specific initial ODE states (e.g. an acclimating tracked
  // state) given the birth environment. No-op for strategies that don't need it.
  void set_initial_states(const environment_type& environment) {
    strategy->set_initial_states(environment, vars);
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
    return strategy->net_mass_production_dt(environment, vars);
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

  // Exact d(growth rate)/d(height) at the current height, delegated to the
  // strategy's AD gradient (#537 A1); returns NA if the strategy provides none,
  // so Node::growth_rate_gradient can fall back to finite differences.
  double growth_rate_gradient_exact(const environment_type& environment) const {
    return strategy->growth_rate_gradient_height_ad(vars.state(HEIGHT_INDEX),
                                                    environment);
  }

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
  basic_internals<S> r_internals() const { return vars; }
  const Control &control() const { return strategy->control; }

private:
  strategy_type_ptr strategy;
  basic_internals<S> vars;
};

template <typename T, typename E> Individual<T,E> make_individual(T s) {
  return Individual<T,E>(make_strategy_ptr(s));
}

} // namespace plant

#endif
