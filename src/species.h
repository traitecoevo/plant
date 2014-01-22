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
#include "interpolator.h"
#include "adaptive_interpolator.h"
#include "state.h"           // get_state, set_state

namespace model {

class SpeciesBase : public ode::OdeTarget {
public:
  typedef std::vector<double> state;
  virtual ~SpeciesBase();
  virtual size_t size() const = 0;
  virtual double height_max() const = 0;
  virtual double leaf_area_above(double height) const = 0;
  virtual std::vector<double> seeds() const = 0;
  virtual void compute_vars_phys(const Environment& environment) = 0;
  virtual void add_seeds(int n) = 0;
  virtual double germination_probability(const Environment&
					 environment) const = 0;
  virtual void clear() = 0;
  virtual void compute_assimilation_fn(const Environment& environment) = 0;
  virtual void rescale_assimilation_fn(const Environment& environment) = 0;
  // Leaf area error estimation
  virtual std::vector<double> leaf_area_error(double scal) const = 0;
  // R-specific wrappers
  virtual std::vector<double> r_height() const = 0;
  virtual void r_set_height(std::vector<double> x) = 0;
  virtual Rcpp::List r_get_plants() const = 0;
  virtual int r_n_individuals() const = 0;
  virtual interpolator::Interpolator r_assimilation_fn() const = 0;
  virtual Strategy r_strategy() const = 0;
  virtual std::vector<double> r_leaf_area() const = 0;
  virtual Rcpp::NumericMatrix r_get_vars_size() const = 0;
  virtual Rcpp::NumericMatrix r_get_vars_phys() const = 0;
  // Still thinking about these...
  virtual state::iterator       get_state(state::iterator it) const = 0;
  virtual state::const_iterator set_state(state::const_iterator it) = 0;
  virtual Rcpp::NumericMatrix r_get_state() const = 0;
  virtual void r_set_state(Rcpp::NumericMatrix x) = 0;
  virtual void r_force_state(Rcpp::NumericMatrix x) = 0;
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
  std::vector<double> seeds() const;
  void compute_vars_phys(const Environment& environment);
  void add_seeds(int n);
  double germination_probability(const Environment& environment) const;
  void clear();
  void compute_assimilation_fn(const Environment& environment);
  void rescale_assimilation_fn(const Environment& environment);

  // * Leaf area error estimation, used by EBT
  std::vector<double> leaf_area_error(double scal) const;

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
  Individual r_seed() const;
  int r_n_individuals() const;
  interpolator::Interpolator r_assimilation_fn() const;
  Strategy r_strategy() const;
  std::vector<double> r_leaf_area() const;
  Rcpp::NumericMatrix r_get_vars_size() const;
  Rcpp::NumericMatrix r_get_vars_phys() const;

  // State
  state::iterator get_state(state::iterator it) const;
  state::const_iterator set_state(state::const_iterator it);
  Rcpp::NumericMatrix r_get_state() const;
  void r_set_state(Rcpp::NumericMatrix x);
  void r_force_state(Rcpp::NumericMatrix x);

private:
  // See Plant for a similar construction
  const Control& control() const {return strategy->get_control();}

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
  for (plants_iterator it = plants.begin();
       it != plants.end(); ++it)
    born += it->offspring();
  return born;
}

// Check to see if individuals have died and remove them from the
// Species.
template <class Individual>
int Species<Individual>::deaths() {
  int died = 0;
  plants_iterator it = plants.begin();
  while (it != plants.end()) {
    if (it->died()) {
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
  for (plants_const_iterator it = plants.begin();
       it != plants.end(); ++it) {
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
std::vector<double> Species<Individual>::seeds() const {
  std::vector<double> ret;
  for (plants_const_iterator it = plants.begin();
       it != plants.end(); ++it) {
    ret.push_back(it->fecundity());
  }
  return ret;
}

// NOTE: We should probably prefer to rescale when this is called
// through the ode stepper.
template <class Individual>
void Species<Individual>::compute_vars_phys(const Environment& environment) {
  if (control().plant_assimilation_approximate_use)
    compute_assimilation_fn(environment);
  for (plants_iterator it = plants.begin();
       it != plants.end(); ++it)
    it->compute_vars_phys(environment);
}
template <>
void Species<CohortTop>::compute_vars_phys(const Environment& environment);

template <class Individual>
void Species<Individual>::add_seeds(int n) {
  for (; n > 0; n--)
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
  std::advance(p, static_cast<int>(util::check_bounds_r(idx, size())));
  return *p;
}

template <class Individual>
Individual Species<Individual>::r_seed() const {
  return seed;
}

template <class Individual>
std::vector<double> Species<Individual>::r_height() const { 
  std::vector<double> ret;
  plants_const_iterator p = plants.begin(); 
  while (p != plants.end())
    ret.push_back((p++)->height());
  return ret;
}

// NOTE: Roll back on error is not possible here at present.
template <class Individual>
void Species<Individual>::r_set_height(std::vector<double> x) {
  util::check_length(x.size(), size());
  if (!util::is_decreasing(x.begin(), x.end()))
    Rcpp::stop("height must be decreasing (ties allowed)");
  std::vector<double>::iterator it = x.begin();
  plants_iterator p = plants.begin();
  while (p != plants.end()) {
    p->set_height(*it++);
    p++;
  }
}

template <class Individual>
Rcpp::List Species<Individual>::r_get_plants() const {
  Rcpp::List ret;
  for (plants_const_iterator it = plants.begin();
       it != plants.end(); ++it)
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
  return static_cast<int>(plants.size());
}
template<> int Species<CohortDiscrete>::r_n_individuals() const;

template <class Individual>
SpeciesBase::state::iterator
Species<Individual>::get_state(SpeciesBase::state::iterator it) const {
  plants_const_iterator first = plants.begin(), last = plants.end();
  while (first != last) {
    it = first->get_state(it);
    ++first;
  }
  return it;
}

template <class Individual>
SpeciesBase::state::const_iterator
Species<Individual>::set_state(SpeciesBase::state::const_iterator it) {
  plants_iterator first = plants.begin(), last = plants.end();
  while (first != last) {
    it = first->set_state(it);
    ++first;
  }
  return it;
}

template <class Individual>
Rcpp::NumericMatrix Species<Individual>::r_get_state() const {
  state tmp(size() * seed.state_size());
  get_state(tmp.begin());

  Rcpp::NumericMatrix ret(static_cast<int>(seed.state_size()),
			  static_cast<int>(size()));
  std::copy(tmp.begin(), tmp.end(), ret.begin());

  return ret;
}

template <class Individual>
void Species<Individual>::r_set_state(Rcpp::NumericMatrix x) {
  util::check_dimensions(static_cast<size_t>(x.nrow()),
			 static_cast<size_t>(x.ncol()),
			 seed.state_size(), size());
  typename Species<Individual>::state tmp(x.begin(), x.end());
  set_state(tmp.begin());
}

// TODO: this does not check that the resulting state is sensible (in
// particular the monotonic check).  That should probably be added.
template <class Individual>
void Species<Individual>::r_force_state(Rcpp::NumericMatrix x) {
  const size_t n = static_cast<size_t>(x.ncol());
  plants.clear();
  for (size_t i = 0; i < n; ++i)
    add_seeds(1);
  r_set_state(x);
}

// Special to Species<CohortTop>
template<> Rcpp::NumericMatrix Species<CohortTop>::r_get_state() const;
template<> void Species<CohortTop>::r_set_state(Rcpp::NumericMatrix x);
template<> void Species<CohortTop>::r_force_state(Rcpp::NumericMatrix x);

// There is some ugliness here; we need to make sure that the smoothed
// unction is not constructed with zero gap between min and max plant
// sizes, so we need to add some epsilon to this.  That is in turn
// duplicated between compute_assimilation_fn and
// rescale_assimilation_fn, for now at least.
template <class Individual>
void Species<Individual>::compute_assimilation_fn(const Environment &environment) {
  const double eps = control().cohort_gradient_eps;
  const double hmin = seed.height() - eps, hmax = height_max() + eps;
  Plant::compute_assimilation_fn(strategy.get(), hmin, hmax,
				     environment);
}

template <class Individual>
void Species<Individual>::rescale_assimilation_fn(const Environment &environment) {
  const double eps = control().cohort_gradient_eps;
  const double hmin = seed.height() - eps, hmax = height_max() + eps;
  Plant::rescale_assimilation_fn(strategy.get(), hmin, hmax,
				 environment);
}

// This doesn't really make any sense for anything other than
// Species<CohortTop>, but the functions are here anyway...
template <class Individual>
std::vector<double> Species<Individual>::leaf_area_error(double scal) const {
  return util::local_error_integration(r_height(), r_leaf_area(), scal);
}

template <class Individual>
interpolator::Interpolator Species<Individual>::r_assimilation_fn() const {
  return strategy->r_assimilation_fn();
}

template <class Individual>
Strategy Species<Individual>::r_strategy() const {
  return *strategy.get();
}

template <class Individual>
std::vector<double> Species<Individual>::r_leaf_area() const {
  std::vector<double> ret;
  plants_const_iterator p = plants.begin();
  while (p != plants.end())
    ret.push_back((p++)->leaf_area());
  return ret;
}

template <class Individual>
Rcpp::NumericMatrix Species<Individual>::r_get_vars_size() const {
  Rcpp::NumericVector tmp = seed.r_get_vars_size();
  Rcpp::NumericMatrix ret(static_cast<int>(tmp.size()),
			  static_cast<int>(size()));

  plants_const_iterator p = plants.begin();
  for (int i = 0; i < ret.ncol(); ++i) {
    ret(Rcpp::_,i) = p->r_get_vars_size();
    ++p;
  }

  Rcpp::CharacterVector names = tmp.names();
  ret.attr("dimnames") = Rcpp::List::create(names, R_NilValue);
  return ret;
}

template <class Individual>
Rcpp::NumericMatrix Species<Individual>::r_get_vars_phys() const {
  Rcpp::NumericVector tmp = seed.r_get_vars_phys();
  Rcpp::NumericMatrix ret(static_cast<int>(tmp.size()),
			  static_cast<int>(size()));

  plants_const_iterator p = plants.begin();
  for (int i = 0; i < ret.ncol(); ++i) {
    ret(Rcpp::_,i) = p->r_get_vars_phys();
    ++p;
  }

  Rcpp::CharacterVector names = tmp.names();
  ret.attr("dimnames") = Rcpp::List::create(names, R_NilValue);
  return ret;
}

template <class Individual>
double Species<Individual>::germination_probability(const Environment& environment) const {
  Individual s(seed);
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
