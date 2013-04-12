// -*-c++-*-
#ifndef TREE_SPECIES_
#define TREE_SPECIES_

#include <list>

#include "ode_target.h"
#include "ode_solver.h"
#include "spline.h"
#include "strategy.h"
#include "cohort_discrete.h" // for a specialisation
#include "util.h" // is_decreasing, check_length

namespace model {

// This class is only used by Patch, and will possibly be where some
// major changes occur once we move to allowing Cohorts.

// It's possible that we should just friend this with Patch.

template <class Individual>
class Species : public ode::OdeTarget {
public:
  Species();
  Species(Strategy *s);

  // * Births and deaths
  int births();
  int deaths();

  // * Lower level functions, used by Patch
  size_t size() const;
  double height_max() const;
  double leaf_area_above(double height) const;
  void compute_vars_phys(spline::Spline *light_environment);
  void add_seeds(int n);
  double germination_probability(spline::Spline *light_environment);
  void clear();

  // * ODE interface
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it, bool &changed);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * R interface
  // Even though there is not an R interface to this class, these
  // functions are used only by other R interface functions.
  std::vector<double> r_get_mass_leaf() const;
  void r_set_mass_leaf(std::vector<double> x);
  Rcpp::List r_get_plants() const;
  int r_n_individuals() const;

private:
  Strategy *strategy;
  Individual seed;
  std::list<Individual> plants;

  typedef typename std::list<Individual>::iterator plants_iterator;
  typedef typename std::list<Individual>::const_iterator plants_const_iterator;
};

// TODO: I'm a bit wary of the impact of seed(NULL) here, especially
// if a copy constructor is triggered, but also just in general.  It
// would be nice if we could just skip this contructor entirely, but
// it is apparently necessary for something in Patch, I think.
template <class Individual>
Species<Individual>::Species() : 
  strategy(NULL),
  seed(strategy) {
}

template <class Individual>
Species<Individual>::Species(Strategy *s) : 
  strategy(s),
  seed(strategy) {
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

template <class Individual>
double Species<Individual>::height_max() const {
  if ( size() == 0 )
    return 0.0;
  else
    return plants.begin()->get_height();
}

template <class Individual>
double Species<Individual>::leaf_area_above(double height) const {
  double tot = 0.0;
  for ( plants_const_iterator it = plants.begin();
	it != plants.end(); it++ )
    tot += it->leaf_area_above(height);
  return tot;
}

template <class Individual>
void Species<Individual>::compute_vars_phys(spline::Spline *light_environment) {
  for ( plants_iterator it = plants.begin();
	it != plants.end(); it++ )
    it->compute_vars_phys(light_environment);
}

template <class Individual>
void Species<Individual>::add_seeds(int n) {
  for ( ; n > 0; n-- )
    plants.push_back(seed);
}

// Declare full specialisation
template<> void Species<CohortDiscrete>::add_seeds(int n);

template <class Individual>
void Species<Individual>::clear() {
  plants.clear();
}

// * ODE interface
template <class Individual>
size_t Species<Individual>::ode_size() const {
  return size() * seed.ode_size();
}

template <class Individual>
ode::iter_const Species<Individual>::ode_values_set(ode::iter_const it, bool &changed) {
  for ( plants_iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->ode_values_set(it, changed);
  return it;
}

template <class Individual>
ode::iter Species<Individual>::ode_values(ode::iter it) const {
  for ( plants_const_iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->ode_values(it);
  return it;
}

template <class Individual>
ode::iter Species<Individual>::ode_rates(ode::iter it) const {
  for ( plants_const_iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->ode_rates(it);
  return it;
}

// * R interface
template <class Individual>
std::vector<double> Species<Individual>::r_get_mass_leaf() const { 
  std::vector<double> ret;
  plants_const_iterator p = plants.begin(); 
  while ( p != plants.end() )
    ret.push_back((p++)->r_get_mass_leaf());
  return ret;
}

// NOTE: Roll back on error is not possible here at present.
template <class Individual>
void Species<Individual>::r_set_mass_leaf(std::vector<double> x) {
  util::check_length(x.size(), size());
  if ( !util::is_decreasing(x.begin(), x.end()) )
    Rf_error("mass_leaf must be decreasing (ties allowed)");
  std::vector<double>::iterator it = x.begin();
  plants_iterator p = plants.begin();
  while ( p != plants.end() ) {
    p->set_mass_leaf(*it++);
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
template<> int Species<CohortDiscrete>::r_n_individuals() const;
template <class Individual>
int Species<Individual>::r_n_individuals() const {
  return (int)plants.size();
}

template <class Individual>
double Species<Individual>::germination_probability(spline::Spline *light_environment) {
  return seed.germination_probability(light_environment);  
}

}

#endif
