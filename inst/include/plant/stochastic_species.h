// -*-c++-*-
#ifndef STOCHASTIC_SPECIES
#define STOCHASTIC_SPECIES

#include <algorithm>
#include <vector>
#include <boost/iterator/filter_iterator.hpp>
#include <plant/util.h>
#include <plant/environment.h>
#include <plant/species_base.h>
#include <plant/stochastic_node.h>
#include <odelia/ode_interface.hpp>

namespace plant {

// This runs the stochastic (finite-population) model: a species is a list of
// discrete individuals that arrive and die as events, rather than a continuous
// size-density distribution discretised into cohorts (that is the deterministic
// Species, in species.h).
//
// The two share their common structure through SpeciesBase (storage + the ODE
// state plumbing + per-element serialisation); see species_base.h. The ODE
// plumbing iterates the *living* individuals here, via the node_begin()/
// node_end() hooks below, which return a filter_iterator over the alive nodes
// (StochasticNode carries its own `alive` flag).
//
// What stays here, because it genuinely differs from the deterministic model:
// competition is a plain sum over individuals (no density weighting, no
// trapezium), and mortality is realised as discrete deaths() that remove
// individuals -- shrinking the ODE system in integer steps -- rather than as a
// continuous mortality state.
//
// Future per-individual tracking (stable ids, birth times) belongs on
// StochasticNode; see the note there and issue #217.
template <typename T, typename E>
class StochasticSpecies
    : public SpeciesBase<StochasticSpecies<T, E>, T, E, StochasticNode<T, E>> {
  typedef SpeciesBase<StochasticSpecies<T, E>, T, E, StochasticNode<T, E>> base_type;
  // Iterators over the living individuals only -- the predicate reads each
  // element's own `alive` flag, so the same odelia free functions the base
  // applies to (always-live) deterministic nodes work here over the subset.
  typedef std::vector<StochasticNode<T, E>> nodes_type;
  typedef boost::filter_iterator<bool (*)(const StochasticNode<T, E>&),
                                  typename nodes_type::iterator> alive_iterator;
  typedef boost::filter_iterator<bool (*)(const StochasticNode<T, E>&),
                                  typename nodes_type::const_iterator> alive_const_iterator;

public:
  typedef T         strategy_type;
  typedef E         environment_type;
  typedef Individual<T,E>  individual_type;
  typedef StochasticNode<T,E> node_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  StochasticSpecies(strategy_type s);

  size_t size() const;
  using base_type::size_individuals;
  void clear();
  void introduce_new_node();
  void introduce_new_node(const E& environment);

  double height_max() const;
  double compute_competition(double height) const;
  void compute_rates(const E& environment);
  double consumption_rate(int i) const;
  std::vector<double> net_reproduction_ratio_by_node() const;

  Rcpp::NumericMatrix r_get_state() const;

  // This is totally new, relative to the deterministic model; this
  // will destructively modify the species by killing individuals.
  size_t deaths();
  double establishment_probability(const E& environment) {
    return new_node.establishment_probability(environment);
  }

  // ODE plumbing is inherited from SpeciesBase; it iterates the living subset
  // through these hooks. (The base resolves them at compile time via CRTP.)
  using base_type::ode_size;
  using base_type::set_ode_state;
  using base_type::ode_state;
  using base_type::ode_rates;
  alive_iterator       node_begin()       { return make_alive(nodes.begin()); }
  alive_iterator       node_end()         { return make_alive(nodes.end()); }
  alive_const_iterator node_begin() const { return make_alive(nodes.begin()); }
  alive_const_iterator node_end()   const { return make_alive(nodes.end()); }

  // * R interface
  std::vector<bool> r_is_alive() const;
  std::vector<double> r_heights() const;
  void r_set_heights(std::vector<double> heights);
  const individual_type& r_new_node() const {return new_node;}
  std::vector<individual_type> r_individuals() const;
  const individual_type& r_individual_at(util::index idx) const {
    return nodes[idx.check_bounds(size_individuals())].individual;
  }

private:
  // Storage lives in the base; bring the members into scope so the many
  // unqualified `nodes` / `strategy` references below resolve.
  using base_type::nodes;
  using base_type::strategy;

  alive_iterator make_alive(typename nodes_type::iterator it) {
    return alive_iterator(&node_type::is_alive, it, nodes.end());
  }
  alive_const_iterator make_alive(typename nodes_type::const_iterator it) const {
    return alive_const_iterator(&node_type::is_alive, it, nodes.end());
  }

  // The template individual for the next introduction; named to mirror
  // Species::new_node (which plays the same role as the next cohort).
  individual_type new_node;
};

template <typename T, typename E>
StochasticSpecies<T,E>::StochasticSpecies(strategy_type s)
  : base_type(s),
    new_node(this->strategy) {
}

template <typename T, typename E>
size_t StochasticSpecies<T,E>::size() const {
  // number of _alive_ plants.
  return std::count_if(nodes.begin(), nodes.end(), &node_type::is_alive);
}

template <typename T, typename E>
void StochasticSpecies<T,E>::clear() {
  nodes.clear();
  // Reset new_node to a blank individual, too.
  new_node = individual_type(strategy);
}

// Note that this does not do establishment probability; suggest that
// this is best to do in the StochasticPatch perhaps?
template <typename T, typename E>
void StochasticSpecies<T,E>::introduce_new_node() {
  nodes.push_back(node_type(new_node));
}

template <typename T, typename E>
void StochasticSpecies<T,E>::introduce_new_node(const E& environment) {
  introduce_new_node();
  // Seed strategy-specific initial states (e.g. TF24's carbohydrate store) from
  // the birth environment before the first rates evaluation, as the
  // deterministic path does in Node::compute_initial_conditions.
  nodes.back().set_initial_states(environment);
  nodes.back().compute_rates(environment);
}


// If a species contains no individuals, we return zero
// (c.f. Species).  Otherwise we return the height of the largest
// individual (always the first in the list).
template <typename T, typename E>
double StochasticSpecies<T,E>::height_max() const {
  for (auto& n : nodes) {
    if (n.alive) {
      return n.height();
    }
  }
  return 0.0;
}

// Because plants are always ordered from largest to smallest, we
// need not continue down the list once the leaf area above a certain
// height is zero, because it will be zero for all plants further down
// the list.
//
// NOTE: In the cases where there is no individuals, we return 0 for
// all heights, as sum(numeric(0)) -> 0
//
// NOTE: A similar early-exit condition to the Plant version is used;
// once the lower bound of the trapezium is zero, we stop including
// individuals.  Working with the boundary plant is tricky here,
// because we might need to include that, too: always in the case of a
// single plant (needed to be the second half of the trapezium) and
// also needed if the last looked at plant was still contributing to
// the integral).
template <typename T, typename E>
double StochasticSpecies<T,E>::compute_competition(double height) const {
  if (size() == 0 || height_max() < height) {
    return 0.0;
  }
  double tot = 0.0;
  for (auto& n : nodes) {
    if (n.alive) {
      if (n.height() > height) {
        tot += n.compute_competition(height);
      } else {
        break;
      }
    }
  }
  return tot;
}

// NOTE: We should probably prefer to rescale when this is called
// through the ode stepper.
template <typename T, typename E>
void StochasticSpecies<T,E>::compute_rates(const E& environment) {
  for (auto& n : nodes) {
    if (n.alive) {
      n.compute_rates(environment);
    }
  }
}

// Uptake of resource `i` summed over the living individuals, each counted once
// (c.f. Species, which integrates the density-weighted rate over the size
// distribution).
template <typename T, typename E>
double StochasticSpecies<T,E>::consumption_rate(int i) const {
  double tot = 0.0;
  for (auto& n : nodes) {
    if (n.alive) {
      tot += n.consumption_rate(i);
    }
  }
  return tot;
}

// TODO(#479): This is going to change...
template <typename T, typename E>
std::vector<double> StochasticSpecies<T,E>::net_reproduction_ratio_by_node() const {
  std::vector<double> ret;
  ret.reserve(size_individuals());
  // I don't think that this is quite right; is it fecundity that we
  // want to track here?  Or do we need to do some more magic to it?
  //
  // basically - I think I need to take the floor here or something?
  //
  // NOTE: dead plants count here!
  for (auto& n : nodes) {
    ret.push_back(n.individual.state(FECUNDITY_INDEX));
  }
  return ret;
}

template <typename T, typename E>
size_t StochasticSpecies<T,E>::deaths() {
  size_t died = 0;
  for (auto& n : nodes) {
    if (n.alive) {
      if (unif_rand() < n.mortality_probability()) {
        n.alive = false;
        died++;
      } else {
        n.reset_mortality();
      }
    }
  }
  return died;
}

template <typename T, typename E>
std::vector<bool> StochasticSpecies<T,E>::r_is_alive() const {
  std::vector<bool> ret;
  ret.reserve(size_individuals());
  for (auto& n : nodes) {
    ret.push_back(n.alive);
  }
  return ret;
}

template <typename T, typename E>
std::vector<Individual<T,E>> StochasticSpecies<T,E>::r_individuals() const {
  std::vector<individual_type> ret;
  ret.reserve(size_individuals());
  for (auto& n : nodes) {
    ret.push_back(n.individual);
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> StochasticSpecies<T,E>::r_heights() const {
  std::vector<double> ret;
  ret.reserve(size());
  // TODO(#479): also simplify r_heights for Species?
  for (auto& n : nodes) {
    if (n.alive) {
      ret.push_back(n.height());
    }
  }
  return ret;
}

template <typename T, typename E>
void StochasticSpecies<T,E>::r_set_heights(std::vector<double> heights) {
  util::check_length(heights.size(), size());
  if (!util::is_decreasing(heights.begin(), heights.end())) {
    util::stop("height must be decreasing (ties allowed)");
  }
  size_t i = 0;
  for (auto& n : nodes) {
    if (n.alive) {
      n.individual.set_state("height", heights[i++]);
    }
  }
}

template <typename T, typename E>
Rcpp::NumericMatrix StochasticSpecies<T,E>::r_get_state() const
{
  size_t ode_size = individual_type::ode_size(), n_individuals = size_individuals();
  size_t aux_size = strategy->aux_size();

  Rcpp::NumericMatrix ret(static_cast<int>(ode_size + aux_size), n_individuals);
  Rcpp::NumericMatrix::iterator it = ret.begin();

  for (size_t i = 0; i < n_individuals; ++i)
  {
    it = this->get_node_state(nodes[i], it);
    it = this->get_node_aux(nodes[i], it);
  }

  // Combine ode_names and aux_names into a single vector for dimnames
  std::vector<std::string> names = individual_type::ode_names();
  std::vector<std::string> aux = strategy->aux_names();
  names.insert(names.end(), aux.begin(), aux.end());

  ret.attr("dimnames") = Rcpp::List::create(names, R_NilValue);
  // Carry per-individual aliveness alongside the state so that a collected
  // snapshot (run_stochastic_collect, via StochasticPatch::state) records which
  // of the (dead-inclusive) columns are still alive. See issue #498.
  ret.attr("is_alive") = Rcpp::wrap(r_is_alive());

  return ret;
}

}

#endif /* STOCHASTIC_SPECIES */
