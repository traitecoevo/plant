// -*-c++-*-
#ifndef TREE_SPECIES_H_
#define TREE_SPECIES_H_

// This one really needs to use some serious Rcpp bits so we're going
// to need to import Rcpp.h I think.  That means being careful about
// where we stick this in tree.h

// TODO: Go through and normalise strategy_type and
// strategy_ptr_type.

// TODO: This file cannot be directly included because it does not
// correctly set up all its dependencies.

#include <vector>
#include <tree2/util.h>

namespace tree2 {

// This is purely for running the deterministic model.

template <typename T>
class Species {
public:
  typedef T cohort_type;
  typedef typename cohort_type::strategy_type strategy_ptr_type;
  typedef typename strategy_ptr_type::element_type strategy_type;
  Species(strategy_type s);

  size_t size() const;
  void clear();
  void add_seed();

  double height_max() const;
  double leaf_area_above(double height) const;
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
  std::vector<double> r_height() const;
  void r_set_height(std::vector<double> height);
  cohort_type r_seed() const {return seed;}
  std::vector<cohort_type> r_plants() const {return plants;}
  // NOTE: I'm in two minds as to whether this is worth exporting.
  cohort_type r_plant_at(util::index idx) const {
    return plants[idx.check_bounds(size())];
  }

  // These are used to determine the degree of cohort refinement.
  std::vector<double> r_leaf_area() const;
  std::vector<double> r_leaf_area_error() const;

private:
  const Control& control() const {return strategy->get_control();}
  strategy_ptr_type strategy;
  cohort_type seed;
  // TODO: this should really say cohorts, rather than plants.
  std::vector<cohort_type> plants;

  typedef typename std::vector<cohort_type>::iterator plants_iterator;
  typedef typename std::vector<cohort_type>::const_iterator plants_const_iterator;
};

template <typename T>
Species<T>::Species(strategy_type s)
  : strategy(make_strategy_ptr(s)),
    seed(strategy) {
}

template <typename T>
size_t Species<T>::size() const {
  return plants.size();
}

template <typename T>
void Species<T>::clear() {
  plants.clear();
  // Reset the seed to a blank seed, too.
  seed = cohort_type(strategy);
}

template <typename T>
void Species<T>::add_seed() {
  plants.push_back(seed);
  // TODO: Should the seed be recomputed here?
}

// If a species contains no individuals, we return the height of a
// seed of the species.  Otherwise we return the height of the largest
// individual (always the first in the list) which will be at least
// tall as a seed.
template <typename T>
double Species<T>::height_max() const {
  return plants.empty() ? seed.height() : plants.front().height();
}

// Because of plants are always ordered from largest to smallest, we
// need not continue down the list once the leaf area above a certain
// height is zero, because it will be zero for all plants further down
// the list.
//
// NOTE: This is simply performing numerical integration, via the
// trapezium rule, of the leaf_area_above with respect to plant
// height.  You'd think that this would be nicer to do in terms of a
// call to an external trapezium integration function, but building
// and discarding the intermediate storage ends up being a nontrivial
// cost.  A more general iterator version might be possible, but with
// the fiddliness around the boundary conditions that won't likely be
// useful.
//
// NOTE: In the cases where there is no individuals, we return 0 for
// all heights.  The integral is not defined, but an empty light
// environment seems appropriate.
//
// NOTE: A similar early-exit condition to the Plant version is used;
// once the lower bound of the trazpeium is zero, we stop including
// individuals.  Working with the boundary cohort is tricky here,
// because we might need to include that, too: always in the case of a
// single cohort (needed to be the second half of the trapezium) and
// also needed if the last looked at plant was still contributing to
// the integral).
template <typename T>
double Species<T>::leaf_area_above(double height) const {
  if (size() == 0 || height_max() < height) {
    return 0.0;
  }
  double tot = 0.0;
  plants_const_iterator it = plants.begin();
  double h1 = it->height(), f_h1 = it->leaf_area_above(height);

  for (++it; it != plants.end(); ++it) {
    const double h0 = it->height(), f_h0 = it->leaf_area_above(height);
    if (!util::is_finite(f_h0)) {
      Rcpp::stop("Detected non-finite contribution");
    }
    tot += (h1 - h0) * (f_h1 + f_h0);
    // Upper point moves for next time:
    h1   = h0;
    f_h1 = f_h0;
    if (h0 < height) {
      break;
    }
  }

  if (size() == 1 || f_h1 > 0) {
    const double h0 = seed.height(), f_h0 = seed.leaf_area_above(height);
    tot += (h1 - h0) * (f_h1 + f_h0);
  }

  return tot / 2;
}

// NOTE: We should probably prefer to rescale when this is called
// through the ode stepper.
template <typename T>
void Species<T>::compute_vars_phys(const Environment& environment) {
  for (auto& p : plants) {
    p.compute_vars_phys(environment);
  }
  seed.compute_initial_conditions(environment);
}

template <typename T>
std::vector<double> Species<T>::seeds() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& p : plants) {
    ret.push_back(p.fecundity());
  }
  return ret;
}

template <typename T>
size_t Species<T>::ode_size() const {
  return size() * T::ode_size();
}

template <typename T>
ode::const_iterator Species<T>::set_ode_state(ode::const_iterator it) {
  return ode::set_ode_state(plants.begin(), plants.end(), it);
}

template <typename T>
ode::iterator Species<T>::ode_state(ode::iterator it) const {
  return ode::ode_state(plants.begin(), plants.end(), it);
}

template <typename T>
ode::iterator Species<T>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(plants.begin(), plants.end(), it);
}


template <typename T>
std::vector<double> Species<T>::r_height() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (plants_const_iterator it = plants.begin();
       it != plants.end(); ++it) {
    ret.push_back(it->height());
  }
  return ret;
}

template <typename T>
void Species<T>::r_set_height(std::vector<double> height) {
  util::check_length(height.size(), size());
  if (!util::is_decreasing(height.begin(), height.end())) {
    util::stop("height must be decreasing (ties allowed)");
  }
  size_t i = 0;
  for (plants_iterator it = plants.begin(); it != plants.end(); ++it, ++i) {
    it->plant.set_height(height[i]);
  }
}

template <typename T>
std::vector<double> Species<T>::r_leaf_area() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& p : plants) {
    ret.push_back(p.leaf_area());
  }
  return ret;
}

template <typename T>
std::vector<double> Species<T>::r_leaf_area_error() const {
  const double scal = 1.0;
  return util::local_error_integration(r_height(), r_leaf_area(), scal);
}

}

#endif
