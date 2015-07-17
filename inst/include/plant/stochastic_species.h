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
// * don't use Cohort<T> for storage
// * support for non-deterministic deaths
// * stochastic birth survival (?)
// * possibly uses a list<T> rather than a vector<T> for storage because
// then we can easily delete without reallocating everything.

// Eventually we might need to support things like stochastic cohorts
// to support:
// * discrete multiple arrivals in a single time
// * tracking birth times
//
// That should be a fairly simple addition here though.
//
// It's possible that this could all be done by some sort of base
// class, or via a composition, but that's not going to happen super
// quickly.
template <typename T>
class StochasticSpecies {
public:
  typedef T plant_type;
  typedef typename T::strategy_type   strategy_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  StochasticSpecies(strategy_type s);

  size_t size() const;
  void clear();
  void add_seed();

  double height_max() const;
  double area_leaf_above(double height) const;
  void compute_vars_phys(const Environment& environment);
  std::vector<double> seeds() const;

  // * ODE interface
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  size_t ode_size() const;
  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;

  // * R interface
  std::vector<double> r_heights() const;
  void r_set_heights(std::vector<double> heights);
  const plant_type& r_seed() const {return seed;}
  std::vector<plant_type> r_plants() const {return plants;}
  const plant_type& r_plant_at(util::index idx) const {
    return plants[idx.check_bounds(size())];
  }

private:
  const Control& control() const {return strategy->get_control();}
  strategy_type_ptr strategy;
  plant_type seed;
  std::vector<plant_type> plants;

  typedef typename std::vector<plant_type>::iterator plants_iterator;
  typedef typename std::vector<plant_type>::const_iterator plants_const_iterator;
};

template <typename T>
StochasticSpecies<T>::StochasticSpecies(strategy_type s)
  : strategy(make_strategy_ptr(s)),
    seed(strategy) {
}

template <typename T>
size_t StochasticSpecies<T>::size() const {
  return plants.size();
}

template <typename T>
void StochasticSpecies<T>::clear() {
  plants.clear();
  // Reset the seed to a blank seed, too.
  seed = plant_type(strategy);
}

// Note that this does not do germination probability; suggest that
// this is best to do in the StochasticPatch perhaps?
template <typename T>
void StochasticSpecies<T>::add_seed() {
  plants.push_back(seed);
}

// If a species contains no individuals, we return the height of a
// seed of the species.  Otherwise we return the height of the largest
// individual (always the first in the list) which will be at least
// tall as a seed.
template <typename T>
double StochasticSpecies<T>::height_max() const {
  return plants.empty() ? seed.height() : plants.front().height();
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
// once the lower bound of the trazpeium is zero, we stop including
// individuals.  Working with the boundary plant is tricky here,
// because we might need to include that, too: always in the case of a
// single plant (needed to be the second half of the trapezium) and
// also needed if the last looked at plant was still contributing to
// the integral).
template <typename T>
double StochasticSpecies<T>::area_leaf_above(double height) const {
  if (size() == 0 || height_max() < height) {
    return 0.0;
  }
  double tot = 0.0;
  for (const auto& it : plants) {
    if (it.height() > height) {
      tot += it.area_leaf_above(height);
    } else {
      break;
    }
  }

  return tot;
}

// NOTE: We should probably prefer to rescale when this is called
// through the ode stepper.
template <typename T>
void StochasticSpecies<T>::compute_vars_phys(const Environment& environment) {
  for (auto& p : plants) {
    p.compute_vars_phys(environment);
  }
  // TODO: this will be worth doing if (and only if) we can store the
  // germination probability.  Otherwise it's probably easiest to do
  // this on demand.
  seed.compute_vars_phys(environment);
}

// TODO: This is going to change...
template <typename T>
std::vector<double> StochasticSpecies<T>::seeds() const {
  std::vector<double> ret;
  ret.reserve(size());
  // I don't think that this is quite right; is it fecundity that we
  // want to track here?  Or do we need to do some more magic to it?
  //
  // basically - I think I need to take the floor here or something?
  for (auto& p : plants) {
    ret.push_back(p.fecundity());
  }
  return ret;
}

template <typename T>
size_t StochasticSpecies<T>::ode_size() const {
  return size() * plant_type::ode_size();
}

template <typename T>
ode::const_iterator StochasticSpecies<T>::set_ode_state(ode::const_iterator it) {
  return ode::set_ode_state(plants.begin(), plants.end(), it);
}

template <typename T>
ode::iterator StochasticSpecies<T>::ode_state(ode::iterator it) const {
  return ode::ode_state(plants.begin(), plants.end(), it);
}

template <typename T>
ode::iterator StochasticSpecies<T>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(plants.begin(), plants.end(), it);
}


template <typename T>
std::vector<double> StochasticSpecies<T>::r_heights() const {
  std::vector<double> ret;
  ret.reserve(size());
  // TODO: also simplify r_heights for Species?
  for (const auto& p : plants) {
    ret.push_back(p.height());
  }
  return ret;
}

template <typename T>
void StochasticSpecies<T>::r_set_heights(std::vector<double> heights) {
  util::check_length(heights.size(), size());
  if (!util::is_decreasing(heights.begin(), heights.end())) {
    util::stop("height must be decreasing (ties allowed)");
  }
  size_t i = 0;
  for (auto& p: plants) {
    p.set_height(heights[i++]);
  }
}

}

#endif
