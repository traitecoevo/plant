// -*-c++-*-
#ifndef PLANT_PLANT_SPECIES_H_
#define PLANT_PLANT_SPECIES_H_

#include <vector>
#include <plant/util.h>
#include <plant/environment.h>
#include <plant/ode_interface.h>
#include <plant/cohort.h>

namespace plant {

// This is purely for running the deterministic model.

template <typename T, typename E>
class Species {
public:
  typedef T         strategy_type;
  typedef E         environment_type;
  typedef Individual<T,E>  individual_type;
  typedef Cohort<T,E> cohort_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Species(strategy_type s);

  size_t size() const;
  void clear();
  void add_seed();

  double height_max() const;
  double compute_competition(double height) const;
  void compute_rates(const environment_type& environment);
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
  const cohort_type& r_seed() const {return seed;}
  std::vector<cohort_type> r_cohorts() const {return cohorts;}
  const cohort_type& r_cohort_at(util::index idx) const {
    return cohorts[idx.check_bounds(size())];
  }

  // Do this with set_ode_state, using an iterator?
  /* double state(int i) const { return vars.state(i); } */
  
  /* double rate(int i) const { return vars.rate(i); } */

  /* void set_state(int i, double v) { */
  /*   vars.set_state(i, v); */
  /* } */


  // These are used to determine the degree of cohort refinement.
  std::vector<double> r_competition_effects() const;
  std::vector<double> r_competition_effects_error(double scal) const;

  // This is just kind of useful
  std::vector<double> r_log_densities() const;

private:
  const Control& control() const {return strategy->get_control();}
  strategy_type_ptr strategy;
  cohort_type seed;
  std::vector<cohort_type> cohorts;

  typedef typename std::vector<cohort_type>::iterator cohorts_iterator;
  typedef typename std::vector<cohort_type>::const_iterator cohorts_const_iterator;
};

template <typename T, typename E>
Species<T,E>::Species(strategy_type s)
  : strategy(make_strategy_ptr(s)),
    seed(strategy) {
}

template <typename T, typename E>
size_t Species<T,E>::size() const {
  return cohorts.size();
}

template <typename T, typename E>
void Species<T,E>::clear() {
  cohorts.clear();
  // Reset the seed to a blank seed, too.
  seed = cohort_type(strategy);
}

template <typename T, typename E>
void Species<T,E>::add_seed() {
  cohorts.push_back(seed);
  // TODO: Should the seed be recomputed here?
}

// If a species contains no individuals, we return the height of a
// seed of the species.  Otherwise we return the height of the largest
// individual (always the first in the list) which will be at least
// tall as a seed.
template <typename T, typename E>
double Species<T,E>::height_max() const {
  return cohorts.empty() ? seed.height() : cohorts.front().height();
}

// Because of cohorts are always ordered from largest to smallest, we
// need not continue down the list once the leaf area above a certain
// height is zero, because it will be zero for all cohorts further down
// the list.
//
// NOTE: This is simply performing numerical integration, via the
// trapezium rule, of the compute_competition with respect to plant
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
template <typename T, typename E>
double Species<T,E>::compute_competition(double height) const {
  if (size() == 0 || height_max() < height) {
    return 0.0;
  }
  double tot = 0.0;
  cohorts_const_iterator it = cohorts.begin();
  double h1 = it->height(), f_h1 = it->compute_competition(height);

  // Loop over cohorts
  for (++it; it != cohorts.end(); ++it) {
    const double h0 = it->height(), f_h0 = it->compute_competition(height);
    if (!util::is_finite(f_h0)) {
      util::stop("Detected non-finite contribution");
    }
    // Integration
    tot += (h1 - h0) * (f_h1 + f_h0);
    // Upper point moves for next time:
    h1   = h0;
    f_h1 = f_h0;
    if (h0 < height) {
      break;
    }
  }

  if (size() == 1 || f_h1 > 0) {
    const double h0 = seed.height(), f_h0 = seed.compute_competition(height);
    tot += (h1 - h0) * (f_h1 + f_h0);
  }

  return tot / 2;
}

// NOTE: We should probably prefer to rescale when this is called
// through the ode stepper.
template <typename T, typename E>
void Species<T,E>::compute_rates(const E& environment) {
  for (auto& c : cohorts) {
    c.compute_rates(environment);
  }
  seed.compute_initial_conditions(environment);
}

template <typename T, typename E>
std::vector<double> Species<T,E>::seeds() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : cohorts) {
    ret.push_back(c.fecundity());
  }
  return ret;
}

template <typename T, typename E>
size_t Species<T,E>::ode_size() const {
  return size() * cohort_type::ode_size();
}

template <typename T, typename E>
ode::const_iterator Species<T,E>::set_ode_state(ode::const_iterator it) {
  return ode::set_ode_state(cohorts.begin(), cohorts.end(), it);
}

template <typename T, typename E>
ode::iterator Species<T,E>::ode_state(ode::iterator it) const {
  return ode::ode_state(cohorts.begin(), cohorts.end(), it);
}

template <typename T, typename E>
ode::iterator Species<T,E>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(cohorts.begin(), cohorts.end(), it);
}
//double sum_aux(int index) {}


template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (cohorts_const_iterator it = cohorts.begin();
       it != cohorts.end(); ++it) {
    ret.push_back(it->height());
  }
  return ret;
}

template <typename T, typename E>
void Species<T,E>::r_set_heights(std::vector<double> heights) {
  util::check_length(heights.size(), size());
  if (!util::is_decreasing(heights.begin(), heights.end())) {
    util::stop("height must be decreasing (ties allowed)");
  }
  size_t i = 0;
  for (cohorts_iterator it = cohorts.begin(); it != cohorts.end(); ++it, ++i) {
    it->plant.set_state("height", heights[i]);
  }
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_competition_effects() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : cohorts) {
    ret.push_back(c.competition_effect());
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_competition_effects_error(double scal) const {
  return util::local_error_integration(r_heights(), r_competition_effects(), scal);
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_densities() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (cohorts_const_iterator it = cohorts.begin();
       it != cohorts.end(); ++it) {
    ret.push_back(it->get_log_density());
  }
  return ret;
}

}

#endif
