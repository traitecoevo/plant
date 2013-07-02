// -*-c++-*-
#ifndef TREE_SPECIES_H_
#define TREE_SPECIES_H_

#include <list>

#include <Rcpp.h>

#include "ode_target.h"
#include "environment.h"
#include "strategy.h"
#include "cohort_discrete.h" // for a specialisation
#include "cohort_top.h"      // for a specialisation
#include "util.h"            // is_decreasing, check_length

namespace model {

class SpeciesBase : public ode::OdeTarget {
public:
  virtual ~SpeciesBase() {};
  virtual size_t size() const = 0;
  virtual double height_max() const = 0;
  virtual double leaf_area_above(double height) const = 0;
  virtual void compute_vars_phys(const Environment& environment) = 0;
  virtual void add_seeds(int n) = 0;
  virtual double germination_probability(const Environment&
					 environment) const = 0;
  virtual void clear() = 0;
  // R-specific wrappers
  virtual std::vector<double> r_height() const = 0;
  virtual void r_set_height(std::vector<double> x) = 0;
  virtual Rcpp::List r_get_plants() const = 0;
  virtual int r_n_individuals() const = 0;
};

template <class Individual>
class Species : public SpeciesBase {
public:
  Species(Strategy  s);
  Species(Strategy *s);

  // * Births and deaths
  int births();
  int deaths();

  // * Lower level functions, used by Patch
  size_t size() const;
  double height_max() const;
  double leaf_area_above(double height) const;
  void compute_vars_phys(const Environment& environment);
  void add_seeds(int n);
  double germination_probability(const Environment& environment) const;
  void clear();

  // * ODE interface
  size_t ode_size() const;
  ode::iterator_const set_ode_values(double time, ode::iterator_const it);
  ode::iterator       ode_values(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it)  const;

  // * R interface
  std::vector<double> r_height() const;
  void r_set_height(std::vector<double> x);
  Rcpp::List r_get_plants() const;
  Individual r_at(size_t idx) const;
  int r_n_individuals() const;
  void r_compute_vars_phys(const Environment& environment);

private:
  Strategy::ptr strategy;
  Individual seed;
  std::list<Individual> plants;

  typedef typename std::list<Individual>::iterator plants_iterator;
  typedef typename std::list<Individual>::const_iterator plants_const_iterator;
};

template <class Individual>
Species<Individual>::Species(Strategy s)
  : strategy(s),
    seed(strategy.get()) {
}

template <class Individual>
Species<Individual>::Species(Strategy *s)
  : strategy(s),
    seed(strategy.get()) {
}

// Compute the number of offspring that will be born by asking all
// individuals how many offspring they will have.
// 
// NOTE: called a second time, this will always return zero, as the
// act of asking plants about how many offspring they have causes them
// to have them.  This may change?
template <class Individual>
int Species<Individual>::births() {
  int born = 0;
  for ( plants_iterator it = plants.begin();
	it != plants.end(); it++ )
    born += it->offspring();
  return born;
}

// Check to see if individuals have died and remove them from the
// Species.
template <class Individual>
int Species<Individual>::deaths() {
  int died = 0;
  plants_iterator it = plants.begin();
  while ( it != plants.end() ) {
    if ( it->died() ) {
      died++;
      it = plants.erase(it); // will advance iterator
    } else {
      it++;
    }
  }

  return died;
}

// * Lower level functions, used by Patch
template <class Individual>
size_t Species<Individual>::size() const {
  return plants.size();
}

// If a species contains no individuals, we return the height of a
// seed of the species.  Otherwise we return the height of the largest
// individual (always the first in the list) which will be at least
// tall as a seed.
template <class Individual>
double Species<Individual>::height_max() const {
  return size() == 0 ? seed.height() : plants.begin()->height();
}

// Because of plants are always ordered from largest to smallest, we
// need not continue down the list once the leaf area above a certain
// height is zero, because it will be zero for all plants further down
// the list.
template <class Individual>
double Species<Individual>::leaf_area_above(double height) const {
  double tot = 0.0;
  for ( plants_const_iterator it = plants.begin();
	it != plants.end(); it++ ) {
    const double a = it->leaf_area_above(height);
    if (a > 0)
      tot += a;
    else
      break;
  }
  return tot;
}
template <>
double Species<CohortTop>::leaf_area_above(double height) const;

template <class Individual>
void Species<Individual>::compute_vars_phys(const Environment& environment) {
  for ( plants_iterator it = plants.begin();
	it != plants.end(); it++ )
    it->compute_vars_phys(environment);
}
template <>
void Species<CohortTop>::compute_vars_phys(const Environment& environment);

template <class Individual>
void Species<Individual>::add_seeds(int n) {
  for ( ; n > 0; n-- )
    plants.push_back(seed);
}

// Declare full specialisation
template <> void Species<CohortDiscrete>::add_seeds(int n);

template <class Individual>
void Species<Individual>::clear() {
  plants.clear();
  // Reset the seed to a blank seed too (relevant for CohortTop,
  // potentially).
  Individual seed_new(strategy.get());
  seed = seed_new;
}

// * ODE interface
template <class Individual>
size_t Species<Individual>::ode_size() const {
  return size() * seed.ode_size();
}

template <class Individual>
ode::iterator_const
Species<Individual>::set_ode_values(double time,
				    ode::iterator_const it) {
  return ode::set_ode_values(plants.begin(), plants.end(), time, it);
}

template <class Individual>
ode::iterator Species<Individual>::ode_values(ode::iterator it) const {
  return ode::ode_values(plants.begin(), plants.end(), it);
}

template <class Individual>
ode::iterator Species<Individual>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(plants.begin(), plants.end(), it);
}

// * R interface
template <class Individual>
Individual Species<Individual>::r_at(size_t idx) const {
  plants_const_iterator p = plants.begin();
  std::advance(p, util::check_bounds_r(idx, size()));
  return *p;
}

template <class Individual>
std::vector<double> Species<Individual>::r_height() const { 
  std::vector<double> ret;
  plants_const_iterator p = plants.begin(); 
  while ( p != plants.end() )
    ret.push_back((p++)->height());
  return ret;
}

// NOTE: Roll back on error is not possible here at present.
template <class Individual>
void Species<Individual>::r_set_height(std::vector<double> x) {
  util::check_length(x.size(), size());
  if ( !util::is_decreasing(x.begin(), x.end()) )
    Rf_error("height must be decreasing (ties allowed)");
  std::vector<double>::iterator it = x.begin();
  plants_iterator p = plants.begin();
  while ( p != plants.end() ) {
    p->set_height(*it++);
    p++;
  }
}

template <class Individual>
Rcpp::List Species<Individual>::r_get_plants() const {
  Rcpp::List ret;
  for ( plants_const_iterator it = plants.begin();
	it != plants.end(); it++ )
    ret.push_back(Rcpp::wrap(*it));
  return ret;
}

// NOTE: This would be slightly easier (avoiding the template
// specialisation) if `Plant` had a method r_n_individuals().
// However, it may not work for other potential Individual classes.
// Also, that converts O(1) -> O(n), which won't be desirable (though
// this is never called for speed).
template <class Individual>
int Species<Individual>::r_n_individuals() const {
  return (int)plants.size();
}
template<> int Species<CohortDiscrete>::r_n_individuals() const;

template <class Individual>
double Species<Individual>::germination_probability(const Environment& environment) const {
  Plant s(seed);
  return s.germination_probability(environment);
}

SEXP species(Rcpp::CppClass individual, Strategy s);

}

// NOTE: I've not chased up why, but I apparently need to use
// RCPP_EXPOSED_CLASS_NOECL here (because these are templated, not
// real), rather than RCPP_EXPOSED_CLASS.
RCPP_EXPOSED_CLASS_NODECL(model::Species<model::Plant>)
RCPP_EXPOSED_CLASS_NODECL(model::Species<model::CohortDiscrete>)
RCPP_EXPOSED_CLASS_NODECL(model::Species<model::CohortTop>)
#endif
