// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_SPECIES_H_
#define PLANT_PLANT_STOCHASTIC_SPECIES_H_

#include <vector>
#include <plant/util.h>
#include <plant/environment.h>
#include <plant/ode_interface.h>

namespace plant {

// This is for running the stochastic model.
//
// It involves a lot of duplication with plant::Species because the
// alternative is some nasty SFINAE / traits work that I don't think
// will help the readability much more than just trying it this way
// first.  We'll see.

// The main difference between this and the deterministic version is that:
// * don't use Cohort<T,E> for storage
// * support for non-deterministic deaths
// * stochastic birth survival (?)

// Eventually we might need to support things like stochastic cohorts
// to support:
// * discrete multiple arrivals in a single time
// * tracking birth times
// * tracking aliveness
//
// That should be a fairly simple addition here though.
//
// It's possible that this could all be done by some sort of base
// class, or via a composition, but that's not going to happen super
// quickly.
template <typename T, typename E>
class StochasticSpecies {
public:
  typedef T         strategy_type;
  typedef E         environment_type;
  typedef Individual<T,E>  individual_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  StochasticSpecies(strategy_type s);

  size_t size() const;
  size_t size_plants() const {return plants.size();}
  void clear();
  void add_seed();
  void add_seed(const E& environment);

  double height_max() const;
  double compute_competition(double height) const;
  void compute_rates(const E& environment);
  std::vector<double> seeds() const;

  // This is totally new, relative to the deterministic model; this
  // will destructively modify the species by removing individuals.
  size_t deaths();
  double establishment_probability(const E& environment) {
    return seed.establishment_probability(environment);
  }

  // * ODE interface
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  size_t ode_size() const;
  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;

  // * R interface
  std::vector<bool> r_is_alive() const {return is_alive;}
  std::vector<double> r_heights() const;
  void r_set_heights(std::vector<double> heights);
  const individual_type& r_seed() const {return seed;}
  std::vector<individual_type> r_plants() const {return plants;}
  const individual_type& r_plant_at(util::index idx) const {
    return plants[idx.check_bounds(size_plants())];
  }

private:
  const Control& control() const {return strategy->get_control();}
  strategy_type_ptr strategy;
  individual_type seed;
  std::vector<individual_type> plants;
  std::vector<bool>       is_alive;
};

template <typename T, typename E>
StochasticSpecies<T,E>::StochasticSpecies(strategy_type s)
  : strategy(make_strategy_ptr(s)),
    seed(strategy) {
}

template <typename T, typename E>
size_t StochasticSpecies<T,E>::size() const {
  // number of _alive_ plants.
  return std::count(is_alive.begin(), is_alive.end(), true);
}

template <typename T, typename E>
void StochasticSpecies<T,E>::clear() {
  plants.clear();
  is_alive.clear();
  // Reset the seed to a blank seed, too.
  seed = individual_type(strategy);
}

// Note that this does not do establishment probability; suggest that
// this is best to do in the StochasticPatch perhaps?
template <typename T, typename E>
void StochasticSpecies<T,E>::add_seed() {
  plants.push_back(seed);
  is_alive.push_back(true);
}

template <typename T, typename E>
void StochasticSpecies<T,E>::add_seed(const E& environment) {
  add_seed();
  plants.back().compute_rates(environment);
}


// If a species contains no individuals, we return zero
// (c.f. Species).  Otherwise we return the height of the largest
// individual (always the first in the list).
template <typename T, typename E>
double StochasticSpecies<T,E>::height_max() const {
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      return plants[i].state(HEIGHT_INDEX);
    }
  }
  return 0.0;
}

// Because of plants are always ordered from largest to smallest, we
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
  // TODO: Here, and elsewhere, consider using a
  // boost::filter_iterator, which is in BH
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      if (plants[i].state(HEIGHT_INDEX) > height) {
        tot += plants[i].compute_competition(height);
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
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      plants[i].compute_rates(environment);
    }
  }
}

// TODO: This is going to change...
template <typename T, typename E>
std::vector<double> StochasticSpecies<T,E>::seeds() const {
  std::vector<double> ret;
  ret.reserve(size());
  // I don't think that this is quite right; is it fecundity that we
  // want to track here?  Or do we need to do some more magic to it?
  //
  // basically - I think I need to take the floor here or something?
  //
  // NOTE: dead plants count here!
  for (auto& p : plants) {
    ret.push_back(p.state(FECUNDITY_INDEX));
  }
  return ret;
}

template <typename T, typename E>
size_t StochasticSpecies<T,E>::deaths() {
  size_t died = 0;
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      if (unif_rand() < plants[i].mortality_probability()) {
        is_alive[i] = false;
        died++;
      } else {
        plants[i].reset_mortality();
      }
    }
  }
  return died;
}


template <typename T, typename E>
size_t StochasticSpecies<T,E>::ode_size() const {
  return size() * individual_type::ode_size();
}

template <typename T, typename E>
ode::const_iterator StochasticSpecies<T,E>::set_ode_state(ode::const_iterator it) {
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      it = plants[i].set_ode_state(it);
    }
  }
  return it;
}

template <typename T, typename E>
ode::iterator StochasticSpecies<T,E>::ode_state(ode::iterator it) const {
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      it = plants[i].ode_state(it);
    }
  }
  return it;
}

template <typename T, typename E>
ode::iterator StochasticSpecies<T,E>::ode_rates(ode::iterator it) const {
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      it = plants[i].ode_rates(it);
    }
  }
  return it;
}


template <typename T, typename E>
std::vector<double> StochasticSpecies<T,E>::r_heights() const {
  std::vector<double> ret;
  ret.reserve(size());
  // TODO: also simplify r_heights for Species?
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      ret.push_back(plants[i].state(HEIGHT_INDEX));
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
  for (size_t i = 0; i < size_plants(); ++i) {
    if (is_alive[i]) {
      plants[i].set_state("height", heights[i]);
    }
  }
}

}

#endif
