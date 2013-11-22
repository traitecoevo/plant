// -*-c++-*-
#ifndef TREE_PATCH_H_
#define TREE_PATCH_H_

#include <vector>
#include <cstddef>

#include "ode_target.h"
#include "adaptive_spline.h"
#include "parameters.h"
#include "species.h"

namespace model {

// The only reason for this is so that the code in interface.cpp does
// not get repeated.
class PatchBase : public ode::OdeTarget {
public:
  virtual ~PatchBase();
  virtual std::vector<int> births() = 0;
  virtual void deaths() = 0;
  virtual size_t size() const = 0;
  virtual double get_time() const = 0;
  virtual void set_time(double x) = 0;
  virtual const Disturbance& get_disturbance_regime() const = 0;
  virtual Rcpp::List r_height() const = 0;
  virtual void r_set_height(Rcpp::List x) = 0;
  virtual Rcpp::List r_get_species() const = 0;
  // r_at not here because it depends on type.
  virtual void r_add_seeds(std::vector<int> seeds) = 0;
  virtual void r_add_seedling(size_t species_index) = 0;
  virtual void r_add_seedlings(std::vector<int> seeds) = 0;
  virtual void reset() = 0;
  virtual std::vector<int> r_n_individuals() const = 0;
  virtual Environment r_environment() const = 0;
  virtual double r_height_max() const = 0;
  virtual double r_canopy_openness(double height) = 0;
  virtual double r_leaf_area_above(double height) const = 0;
  virtual void r_compute_light_environment() = 0;
  virtual void r_compute_vars_phys() = 0;
  virtual std::vector<int> r_germination(std::vector<int> seeds) = 0;
  virtual Parameters r_parameters() const = 0;
  // Still thinking about these
  virtual Rcpp::List r_get_state() const = 0;
  virtual void r_set_state(Rcpp::List x) = 0;
  virtual void r_force_state(Rcpp::List x) = 0;
};

template <class Individual>
class Patch : public PatchBase {
public:
  Patch(Parameters p);
  Patch(Parameters *p);

  // * Lower level functions
  std::vector<int> births();
  void deaths();

  // * Direct manipulation/interrogation
  // Number of species
  size_t size() const;
  void add_seeds(std::vector<int> seeds);
  void add_seedling(size_t species_index);
  void add_seedlings(std::vector<int> seeds);

  // * ODE interface.
  size_t ode_size() const;
  ode::iterator_const set_ode_values(double time, ode::iterator_const it);
  ode::iterator       ode_values(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it)  const;
  double get_time() const;
  void set_time(double x);

  const Disturbance& get_disturbance_regime() const;

  // * R interface.

  // Height (plant size)
  Rcpp::List r_height() const;
  void r_set_height(Rcpp::List x);
  // Access container
  Rcpp::List r_get_species() const;
  const Species<Individual>& at(size_t species_index) const;
  Species<Individual> r_at(size_t species_index) const;
  // Modify container
  void r_add_seeds(std::vector<int> seeds);
  void r_add_seedling(size_t species_index);
  void r_add_seedlings(std::vector<int> seeds);
  void reset();
  // Other interrogation
  std::vector<int> r_n_individuals() const;
  Environment r_environment() const;
  // Wrappers around private methods for use from R
  double r_height_max() const;
  double r_canopy_openness(double height);
  double r_leaf_area_above(double height) const;
  void r_compute_light_environment();
  void r_compute_vars_phys();
  // NOTE: germination is special to the Metacommunity
  std::vector<int> r_germination(std::vector<int> seeds);
  Parameters r_parameters() const;
  // Still thinking about these
  Rcpp::List r_get_state() const;
  void r_set_state(Rcpp::List x);
  void r_force_state(Rcpp::List x);
  
private:
  void initialise();

  // Maximum height for any species in the Patch
  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  double leaf_area_above(double height) const;
  double canopy_openness(double height);

  void compute_light_environment();
  void rescale_light_environment();
  void compute_vars_phys();

  std::vector<int> germination(std::vector<int> seeds);

  void set_state(Rcpp::List x, bool force);

  Parameters::ptr parameters;

  Environment environment;

  std::vector< Species<Individual> > species;

  typedef typename std::vector< Species<Individual> >::iterator 
  species_iterator;
  typedef typename std::vector< Species<Individual> >::const_iterator 
  species_const_iterator;
};

template <class Individual>
Patch<Individual>::Patch(Parameters p)
  : parameters(p),
    environment(*parameters.get()) {
  initialise();
}

template <class Individual>
Patch<Individual>::Patch(Parameters *p)
  : parameters(p),
    environment(*parameters.get()) {
  initialise();
}

template <class Individual>
void Patch<Individual>::deaths() {
  for (species_iterator sp = species.begin();
       sp != species.end(); ++sp)
    sp->deaths();
}

template <class Individual>
std::vector<int> Patch<Individual>::births() {
  std::vector<int> ret(size(), 0);
  std::vector<int>::iterator n = ret.begin();
  for (species_iterator sp = species.begin();
       sp != species.end(); ++sp)
    *n = sp->births();
  return ret;
}

template <class Individual>
void Patch<Individual>::add_seeds(std::vector<int> seeds) {
  seeds = germination(seeds);
  add_seedlings(seeds);
}

template <class Individual>
void Patch<Individual>::add_seedling(size_t species_index) {
  species[species_index].add_seeds(1);
  compute_light_environment();
}

template <class Individual>
void Patch<Individual>::add_seedlings(std::vector<int> seeds) {
  for (size_t i = 0; i < seeds.size(); ++i)
    species[i].add_seeds(seeds[i]);
  compute_light_environment();
}

// There are two components to seed survival:
//   1. survival during dispersal (parameters->Pi_0)
//   2. survival during germination (species-specific, given light env)
// 
// Survival would be distributed as
//   Binom(Binom(seeds_in, p_dispersal), p_germination)
// 
// However, because the dispersal and germination probabilities are
// independent the number that survive is distributed as
//   Binom(seeds_in, p_dispersal * p_germination)
template <class Individual>
std::vector<int> Patch<Individual>::germination(std::vector<int> seeds) {
  const double p_dispersal = parameters->Pi_0;
  for (size_t i = 0; i < seeds.size(); ++i) {
    if (seeds[i] > 0) {
      const double p = p_dispersal *
	species[i].germination_probability(environment);
      seeds[i] = p > 0 ? static_cast<int>(Rf_rbinom(seeds[i], p)) : 0;
    }
  }
  return seeds;
}

// * ODE interface
template <class Individual>
size_t Patch<Individual>::ode_size() const {
  return ode::ode_size(species.begin(), species.end());
}

// NOTE: In theory, recomputing the light environment and
// physiological variables is only necessary if no input variables
// have changed.  This will often be the case on the first call
// (because of the way that derivs() and ode_set_values works, we take
// the values from the model, set them in the ODE solver, then try to
// re-set them in the model.  Importantly, on the first stage of any
// step, nothing has changed because we did get_state -> set_state.
template <class Individual>
ode::iterator_const
Patch<Individual>::set_ode_values(double time,
				  ode::iterator_const it) {
  it = ode::set_ode_values(species.begin(), species.end(), time, it);
  environment.set_time(time);
  if (parameters->control.environment_light_rescale_usually)
    rescale_light_environment();
  else
    compute_light_environment();
  compute_vars_phys();
  return it;
}

template <class Individual>
ode::iterator Patch<Individual>::ode_values(ode::iterator it) const {
  return ode::ode_values(species.begin(), species.end(), it);
}

template <class Individual>
ode::iterator Patch<Individual>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(species.begin(), species.end(), it);
}

template <class Individual>
double Patch<Individual>::get_time() const {
  return environment.get_time();
}
template <class Individual>
void Patch<Individual>::set_time(double x) {
  environment.set_time(x);
}

template <class Individual>
const Disturbance& Patch<Individual>::get_disturbance_regime() const {
  return environment.get_disturbance_regime();
}

// * Private functions
template <class Individual>
void Patch<Individual>::initialise() {
  species.clear(); // (should never be needed?)
  for (size_t i = 0; i < parameters->strategies.size(); ++i) {
    Species<Individual> s(&parameters->strategies[i]);
    species.push_back(s);
  }
  reset();
}

// Number of species
template <class Individual>
size_t Patch<Individual>::size() const {
  return species.size();
}

// Maximum height for any species in the Patch.
//
// Patches with no species have height 0.
// Patches with no individuals have the height of the tallest seedling
// of all species.
template <class Individual>
double Patch<Individual>::height_max() const {
  double ret = 0.0;
  for (species_const_iterator sp = species.begin();
       sp != species.end(); ++sp)
    ret = std::max(ret, sp->height_max());
  return ret;
}

// [eqn 11] Canopy openness at `height`
//
// NOTE: I'd rather that this be a const method (as it is actually
// const) but that conflicts with the definition of DFunctor.
// Probably using Boost and a proper and robust way of binding
// functions would save hassle here.
template <class Individual>
double Patch<Individual>::leaf_area_above(double height) const {
  double tot = 0.0;
  for (species_const_iterator sp = species.begin();
       sp != species.end(); ++sp)
    tot += sp->leaf_area_above(height);
  return tot;
}
template <class Individual>
double Patch<Individual>::canopy_openness(double height) {
  // NOTE: patch_area does not appear in the EBT model formulation.
  return exp(-parameters->c_ext * leaf_area_above(height) /
	     parameters->patch_area);
}

// Create the spline the characterises the light environment.
//
// This needs to be done at least to the height of a seed of the
// tallest species or the height of the tallest individual over all
// species.
template <class Individual>
void Patch<Individual>::compute_light_environment() {
  util::Functor<Patch, &Patch<Individual>::canopy_openness> fun(this);
  environment.compute_light_environment(&fun, height_max());
}

template <class Individual>
void Patch<Individual>::rescale_light_environment() {
  util::Functor<Patch, &Patch<Individual>::canopy_openness> fun(this);
  environment.rescale_light_environment(&fun, height_max());
}

// Given the light environment, "apply" it to every species so that
// physiological variables are updated.
template <class Individual>
void Patch<Individual>::compute_vars_phys() {
  size_t species_index = 0;
  for (species_iterator sp = species.begin();
       sp != species.end(); ++sp) {
    environment.set_seed_rain_index(species_index++);
    sp->compute_vars_phys(environment);
  }
}

// * R interface

// Actually public functions for interrogating & modifying
template <class Individual>
const Species<Individual>& Patch<Individual>::at(size_t species_index) const {
  return species.at(species_index);
}

template <class Individual>
Species<Individual> Patch<Individual>::r_at(size_t species_index) const {
  return species.at(util::check_bounds_r(species_index, size()));
}

template <class Individual>
Rcpp::List Patch<Individual>::r_get_species() const {
  Rcpp::List ret;
  for (species_const_iterator sp = species.begin();
       sp != species.end(); ++sp)
    ret.push_back(Rcpp::wrap(*sp));
  return ret;
}

template <class Individual>
Environment Patch<Individual>::r_environment() const {
  return environment;
}

// In contrast with add_seeds(), we must build the light environment
// in case it has not yet been constructed (when would this be the
// case?).  The post-intruduction light environment is computed via
// add_seeds() -> add_seedlings() -> compute_light_environment().
template <class Individual>
void Patch<Individual>::r_add_seeds(std::vector<int> seeds) {
  util::check_length(seeds.size(), size());
  compute_light_environment();
  add_seeds(seeds);
}

template <class Individual>
void Patch<Individual>::r_add_seedling(size_t species_index) {
  add_seedling(util::check_bounds_r(species_index, size()));
}

template <class Individual>
void Patch<Individual>::r_add_seedlings(std::vector<int> seeds) {
  util::check_length(seeds.size(), size());
  add_seedlings(seeds);
}

// Wrapper functions for testing
template <class Individual>
double Patch<Individual>::r_height_max() const {
  return height_max();
}

template <class Individual>
double Patch<Individual>::r_canopy_openness(double height) {
  return canopy_openness(height);
}

template <class Individual>
double Patch<Individual>::r_leaf_area_above(double height) const {
  return leaf_area_above(height);
}

template <class Individual>
void Patch<Individual>::r_compute_light_environment() {
  compute_light_environment();
}

template <class Individual>
void Patch<Individual>::r_compute_vars_phys() {
  compute_vars_phys();
}

template <class Individual>
std::vector<int> Patch<Individual>::r_germination(std::vector<int> seeds) {
  util::check_length(seeds.size(), size());
  return germination(seeds);
}

template <class Individual>
Parameters Patch<Individual>::r_parameters() const {
  return *parameters.get();
}

template <class Individual>
Rcpp::List Patch<Individual>::r_get_state() const {
  Rcpp::List species_state;
  for (species_const_iterator sp = species.begin();
       sp != species.end(); ++sp)
    species_state.push_back(sp->r_get_state());
  Rcpp::List environment_state = environment.r_get_state();
  return Rcpp::List::create(Rcpp::_["species"]     = species_state,
			    Rcpp::_["environment"] = environment_state);
}

template <class Individual>
void Patch<Individual>::r_set_state(Rcpp::List x) {
  set_state(x, false);
}

template <class Individual>
void Patch<Individual>::r_force_state(Rcpp::List x) {
  set_state(x, true);
}

template <class Individual>
void Patch<Individual>::set_state(Rcpp::List x, bool force) {
  const Rcpp::List species_state = x["species"];
  util::check_length(static_cast<size_t>(species_state.size()), size());
  environment.r_set_state(x["environment"]);
  Rcpp::List::const_iterator it = species_state.begin();
  for (species_iterator sp = species.begin();
       sp != species.end(); ++sp) {
    if (force)
      sp->r_force_state(*it);
    else
      sp->r_set_state(*it);
    ++it;
  }
  compute_light_environment();
  compute_vars_phys();
}

template <class Individual>
Rcpp::List Patch<Individual>::r_height() const {
  Rcpp::List ret;
  for (species_const_iterator sp = species.begin();
       sp != species.end(); ++sp)
    ret.push_back(Rcpp::wrap(sp->r_height()));
  return ret;
}

template <class Individual>
void Patch<Individual>::r_set_height(Rcpp::List x) {
  util::check_length(static_cast<size_t>(x.size()), size());
  for (size_t i = 0; i < size(); ++i)
    species[i].r_set_height(x[static_cast<int>(i)]);

  if (parameters->control.environment_light_rescale_usually)
    rescale_light_environment();
  else
    compute_light_environment();
}

template <class Individual>
std::vector<int> Patch<Individual>::r_n_individuals() const {
  std::vector<int> n;
  for (species_const_iterator sp = species.begin();
       sp != species.end(); ++sp)
    n.push_back(sp->r_n_individuals());
  return n;
}

template <class Individual>
void Patch<Individual>::reset() {
  for (species_iterator sp = species.begin();
       sp != species.end(); ++sp)
    sp->clear();
  environment.clear();
  compute_light_environment();
  compute_vars_phys();
}

SEXP patch(Rcpp::CppClass individual, Parameters p);

}

// NOTE: I've not chased up why, but I apparently need to use
// RCPP_EXPOSED_CLASS_NOECL here (because these are templated, not
// real), rather than RCPP_EXPOSED_CLASS.
RCPP_EXPOSED_CLASS_NODECL(model::Patch<model::Plant>)
RCPP_EXPOSED_CLASS_NODECL(model::Patch<model::CohortDiscrete>)
RCPP_EXPOSED_CLASS_NODECL(model::Patch<model::CohortTop>)

#endif
