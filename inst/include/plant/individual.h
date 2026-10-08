// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_MINIMAL_H_
#define PLANT_PLANT_PLANT_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <odelia/ode_interface.hpp>
#include <vector>
#include <odelia/ode_util.hpp>
#include <plant/census.h>
#include <plant/internals.h>
#include <plant/gradient.h>
#include <plant/util.h>
#include <plant/uniroot.h>
#include <plant/with_slope.h>
#include <utility> // std::pair


namespace plant {

template <typename T, typename E> class Individual {
public:
  using value_type = typename T::value_type;

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
  const value_type& state(std::string name) const {
    return vars.state(strategy->state_index.at(name));
  }
  const value_type& state(int i) const { return vars.state(i); }
  
  // useage:_rate("area_heartwood")
  const value_type& rate(std::string name) const {
    return vars.rate(strategy->state_index.at(name));
  }
  const value_type& rate(int i) const { return vars.rate(i); }

  // useage: set_state("height", 2.0)
  void set_state(std::string name, const value_type& v) {
    int i = strategy->state_index.at(name);
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }
  void set_state(int i, const value_type& v) {
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }

  // aux vars by name and index
  const value_type& aux(std::string name) const {
    return vars.aux(strategy->aux_index.at(name));
  }
  const value_type& aux(int i) const { return vars.aux(i); }

  // set # consumable resources based on env. variables
  void resize_consumption_rates(int i) {
    vars.resize_consumption_rates(i);
  }
  const value_type& consumption_rate(int i) const { return vars.consumption_rate(i); }

  // One census metric read at this individual's own slots. Forwarded to the
  // strategy for the reason every rate function is: the strategy knows which
  // slot holds what and reads it by index.
  // One plant's value of the metric. The reduction over the grid is the
  // species', and the sum over species is the stand's.
  value_type census_value(const census_metric<T>& metric) const {
    return metric.of(*strategy, vars);
  }

  value_type compute_competition(const value_type& z) const {
    return strategy->compute_competition(z, vars);
  }

  // The competition contribution and its vertical derivative, from the one pass
  // the strategy makes. `value` equals compute_competition(z) bit for bit.
  with_slope<value_type>
  compute_competition_and_slope(const value_type& z) const {
    return strategy->compute_competition_and_slope(z, vars);
  }

  // d(leaf area)/d(height), for a caller pulling a leaf-area adjoint back.
  value_type darea_leaf_dheight() const {
    return strategy->darea_leaf_dheight(aux("competition_effect"));
  }

  // The leaf area that adjoint is taken in, for the same caller.
  value_type area_leaf() const { return aux("competition_effect"); }

  // Seed strategy-specific initial ODE states (e.g. an acclimating tracked
  // state) given the birth environment. No-op for strategies that don't need it.
  void set_initial_states(const environment_type& environment) {
    strategy->set_initial_states(environment, vars);
  }

  void compute_rates(const environment_type& environment) {
    if (vars.resource_size() != environment.n_resources()) {
      // handles when Individual hasn't been instantiated in a Patch (ie with an environment)
      vars.resize_consumption_rates(environment.n_resources());
    }
    strategy->compute_rates(environment, vars);
  }
  
  value_type establishment_probability(const environment_type &environment) {
    return strategy->establishment_probability(environment);
  }

  // For a newborn, which sits at birth size: the strategy reads the carbon
  // compute_rates has already left in aux rather than solving the leaf again.
  value_type establishment_probability_of_newborn(const environment_type &environment) {
    return strategy->establishment_probability(environment, vars);
  }

  value_type net_mass_production_dt(const environment_type &environment) {
    // TODO(#483):  maybe reuse intervals? default false
    return strategy->net_mass_production_dt(environment, vars);
  }

  // * ODE interface
  static size_t ode_size() { return strategy_type::state_size(); }
  static std::vector<std::string> ode_names() { return strategy_type::state_names(); }

  // why doesn't strategy->aux_size() also need to be const-qualified? is it because strategy is a pointer?
  size_t aux_size() const { return strategy->aux_size(); }
  std::vector<std::string> aux_names() { return strategy->aux_names(); }

  template <typename It> It set_ode_state(It it) {
    for (size_t i = 0; i < vars.state_size(); i++) {
      vars.states[i] = *it++;
      strategy->update_dependent_aux(i, vars);
    }
    return it;
  }
  template <typename It> It ode_state(It it) const {
    for (size_t i = 0; i < vars.state_size(); i++) {
      util::write_iterator_scalar(it, vars.states[i]);
    }
    return it;
  }
  template <typename It> It ode_rates(It it) const {
    for (size_t i = 0; i < vars.state_size(); i++) {
      util::write_iterator_scalar(it, vars.rates[i]);
    }
    return it;
  }

  template <typename It> It ode_aux(It it) const {
    for (size_t i = 0; i < vars.aux_size(); i++) {
      util::write_iterator_scalar(it, vars.auxs[i]);
    }
    return it;
  }

  // Puts back what ode_aux published, and derives nothing: the slots a state
  // determines are re-derived by set_ode_state, and the rest are quantities a solve
  // produced that this individual is being handed rather than asked to recompute.
  template <typename It> It set_ode_aux(It it) {
    for (size_t i = 0; i < vars.aux_size(); i++) {
      vars.auxs[i] = *it++;
    }
    return it;
  }

  // Its own slots. The strategy is not among them: it is shared, and the species
  // that owns it reaches it once.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, vars);
  }

  // * The recorded block

  // compute_rates as a pure function of its own states, the environment's
  // cohort reads and the strategy's differentiable parameters, in that order.

  // Each section is sized by the class that owns it, so the pack and the
  // adjoint scatter read one layout and cannot drift apart.
  size_t block_input_size(const environment_type& environment) const {
    return strategy_type::state_size() + environment.n_cohort_reads() +
           strategy->ad_parameters().size();
  }

  // The strategy rates, then the density transport term, then one consumption
  // rate per resource.
  size_t block_output_size(const environment_type& environment) const {
    return strategy_type::state_size() + 1 + environment.n_resources();
  }

  template <typename It>
  It block_inputs(It it, const environment_type& environment) const {
    it = ode_state(it);
    it = environment.cohort_reads(it);
    for (const value_type* p : strategy->ad_parameters()) {
      util::write_iterator_scalar(it, *p);
    }
    return it;
  }

  // The states go in through set_state, never into vars.states: only
  // update_dependent_aux writes competition_effect and height_inverse.

  // They are applied last, after the parameters, because area_leaf(height)
  // reads lma and would otherwise be derived at the previous block's.
  template <typename It>
  It set_block_inputs(It it, environment_type& environment) {
    std::vector<value_type> state(vars.state_size());
    for (size_t i = 0; i < vars.state_size(); ++i) {
      state[i] = *it++;
    }
    it = environment.set_cohort_reads(it);
    for (value_type* p : strategy->ad_parameters()) {
      *p = *it++;
    }
    for (size_t i = 0; i < vars.state_size(); ++i) {
      set_state(static_cast<int>(i), state[i]);
    }
    return it;
  }

  // The transport term is read after the rates and before the consumption
  // rates: it runs a second rate evaluation on a copy, and the copy shares this
  // strategy, so the rates already read must be off the strategy first.
  template <typename It>
  It block_outputs(It it, const environment_type& environment) const {
    it = ode_rates(it);
    util::write_iterator_scalar(it, log_density_rate(environment));
    for (size_t i = 0; i < vars.resource_size(); ++i) {
      util::write_iterator_scalar(it, vars.consumption_rates[i]);
    }
    return it;
  }

  // The size-density equation's transport term, from the rates just computed.
  //
  // A density in birth date changes only by mortality: nothing moves an
  // individual along the birth-date axis, so there is no compression term and
  // no displaced-height probe to pay for. A density in height additionally
  // compresses as the spacing between neighbouring sizes changes. The branch is
  // here rather than in Node because this rate is one of the block's outputs,
  // and a block that emitted the height-coordinate term on the birth-date
  // coordinate would be differentiating a term that is identically zero.
  value_type log_density_rate(const environment_type& environment) const {
    if (control().node_density_in_birth_date) {
      return -rate(MORTALITY_INDEX);
    }
    return -growth_rate_gradient(environment) - rate(MORTALITY_INDEX);
  }

  // d(growth rate)/d(height), as a sub-grid probe: one further rate evaluation
  // at a displaced height, differenced against the rate already computed.
  // Differencing the growth rate needs a mutable Individual to perturb height
  // on, and it must not disturb this one's computed state and rates, so perturb
  // a copy. Both evaluations carry the scalar, so a recording of this holds the
  // quotient and the reverse pass forms no seed for the growth rate by hand.
  value_type growth_rate_gradient(const environment_type& environment) const {
    Individual probe = *this;
    // The lambda carries value_type in and out. Written as double it would
    // still compile, taking the value of an active growth rate, and the
    // transport term would then be built from a derivative of exactly zero with
    // nothing raised.
    auto fun = [&](const value_type& h) -> value_type {
      return probe.growth_rate_given_height(h, environment);
    };
    const Control& ctrl = control();
    const double eps = ctrl.node_gradient_eps;
    if (ctrl.node_gradient_richardson) {
      return util::gradient_richardson(fun, state(HEIGHT_INDEX), eps,
                                       ctrl.node_gradient_richardson_depth);
    }
    return util::gradient_fd(fun, state(HEIGHT_INDEX), eps,
                             rate(HEIGHT_INDEX), ctrl.node_gradient_direction);
  }

  // Single individual methods

  // Used in the stochastic model:
  value_type mortality_probability() const { return 1 - exp(-state(MORTALITY_INDEX)); }
  
  void reset_mortality() { set_state("mortality", 0.0); }

  // Carries value_type in both directions. Declaring either side double would
  // still compile, because the value would be taken through a return type
  // written as double, and Node::growth_rate_gradient differences this function:
  // the transport term would then read a derivative of exactly zero with nothing
  // raised, and log_density_dt is built from it.
  value_type growth_rate_given_height(const value_type& height,
                                      const environment_type& environment) {
    // Called repeatedly from the finite-difference gradient (Node::
    // growth_rate_gradient), so address height by integer slot rather than the
    // "height" string-map lookup (see #466).
    set_state(HEIGHT_INDEX, height);
    compute_rates(environment);
    return rate(HEIGHT_INDEX);
  }

  // Reached only from R, which takes a double, and the point is located by
  // iterating on the residual rather than by a declared derivative, so the
  // residual is read at its value here.
  double resource_compensation_point() {
    environment_type env = environment_type();

    auto target = [&] (double x) mutable -> double {
      env.set_fixed_environment(x, 100);
      compute_rates(env);
      return odelia::util::to_passive(net_mass_production_dt(env));
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
  // R takes a double, so the store is read out at its values.
  Internals<double> r_internals() const {
    Internals<double> out(vars.state_size(), vars.aux_size(), vars.resource_size());
    for (size_t i = 0; i < vars.state_size(); ++i) {
      out.states[i] = odelia::util::to_passive(vars.states[i]);
      out.rates[i]  = odelia::util::to_passive(vars.rates[i]);
    }
    for (size_t i = 0; i < vars.aux_size(); ++i) {
      out.auxs[i] = odelia::util::to_passive(vars.auxs[i]);
    }
    for (size_t i = 0; i < vars.resource_size(); ++i) {
      out.consumption_rates[i] =
        odelia::util::to_passive(vars.consumption_rates[i]);
    }
    return out;
  }
  const Control &control() const { return strategy->control; }

private:
  strategy_type_ptr strategy;
  Internals<value_type> vars;
};

template <typename T, typename E> Individual<T,E> make_individual(T s) {
  return Individual<T,E>(make_strategy_ptr(s));
}

} // namespace plant

#endif
