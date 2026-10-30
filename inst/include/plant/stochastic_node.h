// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_NODE_H_
#define PLANT_PLANT_STOCHASTIC_NODE_H_

#include <plant/individual.h>
#include <odelia/ode_interface.hpp>
#include <plant/with_slope.h>

namespace plant {

// The storage element of StochasticSpecies: a single finite individual plus the
// bookkeeping the stochastic model needs but the deterministic Node does not.
//
// This is the counterpart to Node<T,E> for the finite-population solver. Where
// Node augments an Individual with the continuous size-density machinery
// (log_density, survival-weighted offspring, birth bookkeeping), StochasticNode
// augments it with a discrete-population concern: whether the individual is
// still alive. Folding `alive` into the element (rather than a vector<bool> held
// parallel to the individuals, as before) lets a predicate decide liveness from
// the element alone, so the same odelia ODE free functions that serve Species
// can iterate the *living* subset here -- see StochasticSpecies.
//
// The ODE state of a StochasticNode is exactly its Individual's state: unlike
// Node there are no extra density/offspring equations. The methods below simply
// forward to the Individual so the element satisfies odelia's element interface
// (ode_size / set_ode_state / ode_state / ode_rates).
//
// This is also the natural home for the per-individual tracking the stochastic
// runner has long wanted (issue #217, stochastic_patch_runner.h): a stable id
// and a birth time, so a death cannot scramble which trajectory is which. Those
// are not added yet, but they belong here when they are.
template <typename T, typename E>
class StochasticNode {
public:
  using value_type = typename T::value_type;

  typedef Individual<T, E> individual_type;

  explicit StochasticNode(individual_type individual_)
    : individual(individual_), alive(true) {}

  // --- forwards used by StochasticSpecies / StochasticPatch ---------------
  value_type height() const { return individual.state(HEIGHT_INDEX); }
  value_type compute_competition(const value_type& z) const {
    return individual.compute_competition(z);
  }
  with_slope<value_type>
  compute_competition_and_slope(const value_type& z) const {
    return individual.compute_competition_and_slope(z);
  }
  void compute_rates(const E& environment) {
    individual.compute_rates(environment);
  }
  void set_initial_states(const E& environment) {
    individual.set_initial_states(environment);
  }
  value_type mortality_probability() const {
    return individual.mortality_probability();
  }
  void reset_mortality() { individual.reset_mortality(); }
  // One individual's uptake of resource `i`. The deterministic Node scales this
  // by cohort density; here the individual is the unit.
  double consumption_rate(int i) const {
    return individual.consumption_rate(i);
  }

  // --- odelia element ODE interface (state only; no density/offspring) ----
  static size_t ode_size() { return individual_type::ode_size(); }
  size_t aux_size() const { return individual.aux_size(); }

  template <typename It> It set_ode_state(It it) {
    return individual.set_ode_state(it);
  }
  template <typename It> It ode_state(It it) const {
    return individual.ode_state(it);
  }
  template <typename It> It ode_rates(It it) const {
    return individual.ode_rates(it);
  }
  template <typename It> It ode_aux(It it) const {
    return individual.ode_aux(it);
  }

  // Predicate for filtering the living subset (see StochasticSpecies). A free
  // function on the element keeps the filter usable by boost::filter_iterator.
  static bool is_alive(const StochasticNode<T, E>& n) { return n.alive; }

  individual_type individual;
  bool alive;
};

} // namespace plant

#endif
