// -*-c++-*-
#ifndef STOCHASTIC_SPECIES
#define STOCHASTIC_SPECIES

#include <algorithm>
#include <vector>
#include <boost/iterator/filter_iterator.hpp>
#include <plant/util.h>
#include <plant/models/environment.h>
#include <plant/community/species_base.h>
#include <plant/community/stochastic_node.h>
#include <odelia/ode_interface.hpp>
#include <plant/numerics/with_slope.h>

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
  using value_type = typename T::value_type;

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

  value_type height_max() const;
  value_type compute_competition(double height) const;
  // The same sum and its vertical derivative, from one pass over the living
  // individuals. `value` equals compute_competition(height) bit for bit.
  with_slope<value_type>
  compute_competition_and_slope(double height) const;
  void compute_rates(const E& environment);
  double consumption_rate(int i) const;
  std::vector<double> net_reproduction_ratio_by_node() const;

  Rcpp::NumericMatrix r_get_state() const;

  // This is totally new, relative to the deterministic model; this
  // will destructively modify the species by killing individuals.
  size_t deaths();
  // An event's removal: each living individual is removed with probability
  // pr(height), the stochastic counterpart of the deterministic patch scaling a
  // cohort's density by 1 - pr. One with pr(height) = 0 is not drawn for:
  // an event takes no random number for an individual it cannot reach.
  // Returns {individuals drawn for, individuals removed}.
  template <typename Pr>
  std::pair<size_t, size_t> remove(Pr pr);
  // Arrival times in (from, to) on `area`: a Poisson process whose rate is
  // the strategy's birth rate per unit area, the same curve the deterministic
  // solver's inflow reads. Drawn exactly, by thinning: on each piece between
  // the rate's knots (one piece when it is constant), candidates arrive at the
  // piece's largest rate, by exponential gaps, and each is kept with
  // probability rate / that. Between its knots a varying rate is a cubic, so
  // four values fix it, and its largest is at an end or at a turning point of
  // that cubic; the spline being monotone there is not assumed. Strictly
  // ascending, which a schedule needs, since one species is introduced once
  // at a time.
  std::vector<double> draw_arrivals(double from, double to, double area) const;
  // More arrivals than this from one draw is refused rather than allocated.
  static constexpr double max_arrivals = 1e7;
  value_type establishment_probability(const E& environment) {
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
  // Initialise strategy-specific initial states (e.g. TF24's carbohydrate store) from
  // the birth environment before the first rates evaluation, as the
  // deterministic path does in Node::compute_initial_conditions.
  nodes.back().set_initial_states(environment);
  nodes.back().compute_rates(environment);
}


// The height of the tallest living individual, or zero if there is none (c.f.
// Species). Individuals are kept in order of arrival, which is not an order of
// height in a model whose plants can overtake older ones, so every one is read.
template <typename T, typename E>
typename StochasticSpecies<T,E>::value_type
StochasticSpecies<T,E>::height_max() const {
  value_type ret = 0.0;
  for (auto& n : nodes) {
    if (n.alive && n.height() > ret) {
      ret = n.height();
    }
  }
  return ret;
}

// Leaf area above `height`: the sum over every living individual taller than
// it, in order of arrival. Zero for an empty species.
template <typename T, typename E>
typename StochasticSpecies<T,E>::value_type
StochasticSpecies<T,E>::compute_competition(double height) const {
  value_type tot = 0.0;
  for (auto& n : nodes) {
    if (n.alive && n.height() > height) {
      tot += n.compute_competition(height);
    }
  }
  return tot;
}

template <typename T, typename E>
with_slope<typename StochasticSpecies<T,E>::value_type>
StochasticSpecies<T,E>::compute_competition_and_slope(double height) const {
  with_slope<value_type> sum{0.0, 0.0};
  for (auto& n : nodes) {
    if (n.alive && n.height() > height) {
      const with_slope<value_type> fs =
        n.compute_competition_and_slope(height);
      sum.value += fs.value;
      sum.slope += fs.slope;
    }
  }
  return sum;
}

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

// Each individual's fecundity state, the dead included at the value they died
// with (only the living are integrated). This is
// not yet a net reproduction ratio for the stochastic model: what that should
// be, and whether an individual that has died still counts, is open.
template <typename T, typename E>
std::vector<double> StochasticSpecies<T,E>::net_reproduction_ratio_by_node() const {
  std::vector<double> ret;
  ret.reserve(size_individuals());
  for (auto& n : nodes) {
    ret.push_back(
      odelia::util::to_passive(n.individual.state(FECUNDITY_INDEX)));
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
std::vector<double>
StochasticSpecies<T,E>::draw_arrivals(double from, double to,
                                      double area) const {
  const auto& s = *this->strategy;
  const bool varies = s.is_variable_birth_rate;
  auto rate = [&](double t) {
    return odelia::util::to_passive(s.birth_rate_at(t));
  };
  // The largest rate on [a, b]: the cubic through the rate at a, a + h/3,
  // a + 2h/3 and b, written in s = 3 (t - a) / h by forward differences, is
  // y0 + A s + B s^2 + C s^3; its turning points are where A + 2 B s + 3 C s^2
  // vanishes, and the rate is read there rather than off the fit.
  auto largest = [&](double a, double b) {
    const double h = b - a;
    const double y0 = rate(a), y1 = rate(a + h / 3), y2 = rate(a + 2 * h / 3),
                 y3 = rate(b);
    double top = std::max(std::max(y0, y1), std::max(y2, y3));
    if (!varies) {
      return top;
    }
    const double d1 = y1 - y0, d2 = y2 - 2 * y1 + y0,
                 d3 = y3 - 3 * y2 + 3 * y1 - y0;
    const double A = d1 - d2 / 2 + d3 / 3, B = d2 / 2 - d3 / 2, C = d3 / 6;
    std::vector<double> roots;
    if (C != 0.0) {
      const double disc = 4 * B * B - 12 * A * C;
      if (disc >= 0.0) {
        roots.push_back((-2 * B + std::sqrt(disc)) / (6 * C));
        roots.push_back((-2 * B - std::sqrt(disc)) / (6 * C));
      }
    } else if (B != 0.0) {
      roots.push_back(-A / (2 * B));
    }
    for (const double r : roots) {
      if (r > 0.0 && r < 3.0) {
        top = std::max(top, rate(a + r * h / 3));
      }
    }
    return top;
  };

  std::vector<double> edges = {from};
  if (varies) {
    for (const double k : s.extrinsic_drivers.knots("birth_rate", from, to)) {
      edges.push_back(k);
    }
  }
  edges.push_back(to);
  std::vector<double> ret;
  for (size_t i = 1; i < edges.size(); ++i) {
    const double a = edges[i - 1], b = edges[i];
    if (!(b > a)) {
      continue;
    }
    const double bound = largest(a, b);
    if (!std::isfinite(bound)) {
      util::stop("The birth rate is not finite, so no arrivals can be drawn");
    }
    if (!(bound > 0.0)) {
      continue;
    }
    const double per_year = area * bound;
    double t = a;
    for (;;) {
      t += exp_rand() / per_year;
      if (!(t < b)) {
        break;
      }
      if (varies) {
        const double r = rate(t);
        // The bound is the cubic's largest, so this is rounding at most; more
        // would make the draw an approximation, silently.
        if (r > bound * (1.0 + 1e-9)) {
          util::stop("The birth rate exceeds its bound between knots, so "
                     "arrivals cannot be drawn exactly");
        }
        if (!(unif_rand() * bound < r)) {
          continue;
        }
      }
      // A gap below the time's rounding would repeat it.
      if (!ret.empty() && !(t > ret.back())) {
        continue;
      }
      ret.push_back(t);
      if (static_cast<double>(ret.size()) > max_arrivals) {
        util::stop("More than " + util::to_string(max_arrivals) +
                   " arrivals in one draw; lower patch_area or the birth rate");
      }
    }
  }
  return ret;
}

template <typename T, typename E>
template <typename Pr>
std::pair<size_t, size_t> StochasticSpecies<T,E>::remove(Pr pr) {
  size_t drawn = 0, removed = 0;
  for (auto& n : nodes) {
    if (!n.alive) {
      continue;
    }
    const double p = pr(odelia::util::to_passive(n.height()));
    if (!(p > 0.0)) {
      continue;
    }
    ++drawn;
    if (unif_rand() < p) {
      n.alive = false;
      ++removed;
    }
  }
  return {drawn, removed};
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
  for (auto& n : nodes) {
    if (n.alive) {
      ret.push_back(odelia::util::to_passive(n.height()));
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
  // snapshot (a stochastic run's record = "steps", via StochasticPatch::state)
  // records which of the (dead-inclusive) columns are still alive. See issue
  // #498.
  ret.attr("is_alive") = Rcpp::wrap(r_is_alive());

  return ret;
}

}

#endif /* STOCHASTIC_SPECIES */
