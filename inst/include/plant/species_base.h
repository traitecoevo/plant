// -*-c++-*-
#ifndef PLANT_PLANT_SPECIES_BASE_H_
#define PLANT_PLANT_SPECIES_BASE_H_

#include <vector>
#include <plant/individual.h>
#include <odelia/ode_interface.hpp>

namespace plant {

// Common structure shared by the two species containers:
//   * Species<T,E>            -- the deterministic size-density distribution,
//                               discretised into cohorts (Node elements); and
//   * StochasticSpecies<T,E>  -- the finite population of discrete individuals
//                               (StochasticNode elements).
//
// Both own a strategy and a vector of an element type, and both present the
// same ODE-object interface to the patch above them. This CRTP base owns that
// common machinery -- the strategy/element storage, the ODE state plumbing, and
// the per-element (de)serialisation helpers -- so it is written once rather than
// duplicated. The element type is a template parameter (Node or StochasticNode);
// both expose the same `individual` member and odelia element ODE interface.
//
// What is deliberately NOT here is the model logic that genuinely differs and
// sits on (or near) the hot path: the competition contribution (a density-
// weighted trapezium for the deterministic model, a plain sum over living
// individuals for the stochastic one) and the rate computation (survival-
// weighted, with boundary-node initial conditions, vs a bare per-individual
// step plus discrete deaths). Those stay in the derived classes, defined
// directly on the concrete element type, so the codebase's hot-path performance
// choices (agents.md s12) are unaffected by this sharing -- there is no runtime
// dispatch, and the divergent methods inline exactly as before.
//
// The ODE plumbing iterates a range the derived class supplies via
// node_begin()/node_end(): all nodes for the deterministic Species, the living
// subset (a filter_iterator) for StochasticSpecies. Those hooks are resolved at
// compile time through the CRTP downcast.
template <typename Derived, typename T, typename E, typename Element>
class SpeciesBase {
public:
  typedef T                           strategy_type;
  typedef E                           environment_type;
  typedef Individual<T, E>            individual_type;
  typedef Element                     node_type;
  typedef typename strategy_type::ptr strategy_type_ptr;

  using value_type = typename T::value_type;

  size_t size_individuals() const { return nodes.size(); }

  // * ODE interface -- delegates to the shared odelia free functions over the
  // derived class's chosen element range (all nodes, or the living subset).
  size_t ode_size() const {
    return odelia::ode::ode_size(d().node_begin(), d().node_end());
  }
  template <typename It> It set_ode_state(It it) {
    return odelia::ode::set_ode_state(d().node_begin(), d().node_end(), it);
  }
  template <typename It> It ode_state(It it) const {
    return odelia::ode::ode_state(d().node_begin(), d().node_end(), it);
  }
  template <typename It> It ode_rates(It it) const {
    return odelia::ode::ode_rates(d().node_begin(), d().node_end(), it);
  }
  template <typename It> It ode_aux(It it) const {
    return odelia::ode::ode_aux(d().node_begin(), d().node_end(), it);
  }
  template <typename It> It set_ode_aux(It it) {
    return odelia::ode::set_ode_aux(d().node_begin(), d().node_end(), it);
  }

  // Serialise one element's ODE state / aux into an R matrix column. Generic
  // over the element: Node carries the extra density/offspring equations, a
  // StochasticNode forwards to its Individual, so each writes the right rows.
  Rcpp::NumericMatrix::iterator
  get_node_state(const node_type& node, Rcpp::NumericMatrix::iterator it) const {
    std::vector<double> tmp = odelia::ode::r_ode_state(node);
    return std::copy(tmp.begin(), tmp.end(), it);
  }
  Rcpp::NumericMatrix::iterator
  get_node_aux(const node_type& node, Rcpp::NumericMatrix::iterator it) const {
    std::vector<double> tmp = odelia::ode::r_ode_aux(node);
    return std::copy(tmp.begin(), tmp.end(), it);
  }

protected:
  explicit SpeciesBase(strategy_type s) : strategy(make_strategy_ptr(s)) {}
  // An already-prepared strategy, taken as is. prepare_strategy() is refused at
  // an active scalar, so a species there is reachable only through this.
  explicit SpeciesBase(strategy_type_ptr s) : strategy(s) {}

  const Control& control() const { return strategy->control; }

  strategy_type_ptr strategy;
  std::vector<node_type> nodes;

private:
  Derived&       d()       { return static_cast<Derived&>(*this); }
  const Derived& d() const { return static_cast<const Derived&>(*this); }
};

} // namespace plant

#endif
