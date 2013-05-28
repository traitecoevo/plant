// -*-c++-*-
#ifndef TREE_PATCH_H_
#define TREE_PATCH_H_

#include <vector>

#include "ode_target.h"
#include "adaptive_spline.h"
#include "parameters.h"
#include "species.h"

namespace model {

// The only reason for this is so that the code in interface.cpp does
// not get repeated.
class PatchBase : public ode::OdeTarget {
public:
  virtual ~PatchBase() {};
  virtual void r_step() = 0;
  virtual void step_deterministic() = 0;
  virtual void r_step_stochastic() = 0;
  virtual size_t size() const = 0;
  virtual void r_add_seeds(std::vector<int> seeds) = 0;
  virtual void r_add_seedlings(std::vector<int> seeds) = 0;
  virtual Rcpp::List r_height() const = 0;
  virtual void r_set_height(Rcpp::List x) = 0;
  virtual Rcpp::List r_get_species() const = 0;
  virtual void clear() = 0;
  virtual double r_age() const = 0;
  virtual std::vector<int> r_n_individuals() const = 0;
  virtual double r_height_max() const = 0;
  virtual Environment r_environment() const = 0;
  virtual double r_canopy_openness(double height) = 0;
  virtual void r_compute_light_environment() = 0;
  virtual void r_compute_vars_phys() = 0;
  virtual std::vector<int> r_germination(std::vector<int> seeds) = 0;
};

template <class Individual>
class Patch : public PatchBase {
public:
  Patch(Parameters p);
  Patch(Parameters *p);

  // * Methods used from Metacommunity:
  // Advance the system through one complete time step.
  void step();
  // Advance the system through one time step deterministically
  // (plant growth, physiological accounting)
  void step_deterministic();
  // Advance the system through the stochastic life cycle stages
  // (producing seeds and dying).
  void step_stochastic();

  // * Lower level functions
  std::vector<int> births();
  void deaths();

  // * Direct manipulation/interrogation
  // Number of species
  size_t size() const;
  void add_seeds(std::vector<int> seeds);
  void add_seedlings(std::vector<int> seeds);

  // * ODE interface.
  void derivs(double time, ode::iter_const y, ode::iter dydt);
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * R interface.

  // Height (plant size)
  Rcpp::List r_height() const;
  void r_set_height(Rcpp::List x);
  // Access container
  Rcpp::List r_get_species() const;
  Species<Individual> r_at(size_t idx) const;
  // Modify container
  void r_add_seeds(std::vector<int> seeds);
  void r_add_seedlings(std::vector<int> seeds);
  void clear();
  // Other interrogation
  double r_age() const {return environment.get_age();}
  std::vector<int> r_n_individuals() const;
  Environment r_environment() const;
  void r_step();
  void r_step_stochastic();
  double r_height_max() const { return height_max(); }
  double r_canopy_openness(double height) {return canopy_openness(height);}
  // Wrappers for testing
  void r_compute_light_environment() {compute_light_environment();}
  void r_compute_vars_phys() {compute_vars_phys();}
  std::vector<int> r_germination(std::vector<int> seeds);
  
private:
  void initialise();

  // Maximum height for any species in the Patch
  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  double canopy_openness(double height);

  void compute_light_environment();
  void compute_vars_phys();

  std::vector<int> germination(std::vector<int> seeds);

  Parameters::ptr parameters;

  Environment environment;

  std::vector< Species<Individual> > species;
  ode::Solver<Patch> ode_solver;

  typedef typename std::vector< Species<Individual> >::iterator 
  species_iterator;
  typedef typename std::vector< Species<Individual> >::const_iterator 
  species_const_iterator;
};

template <class Individual>
Patch<Individual>::Patch(Parameters p)
  : parameters(p),
    environment(*parameters.get()),
    ode_solver(this) {
  initialise();
}

template <class Individual>
Patch<Individual>::Patch(Parameters *p)
  : parameters(p),
    environment(*parameters.get()),
    ode_solver(this) {
  initialise();
}

template <class Individual>
void Patch<Individual>::step() {
  step_deterministic();
  step_stochastic();
}

template <class Individual>
void Patch<Individual>::step_deterministic() {
  std::vector<double> y(ode_size());
  ode_values(y.begin());
  ode_solver.set_state(y, environment.get_age());
  ode_solver.step();
  environment.set_age(ode_solver.get_time());
}

// Note that this will change when we are working with a
// metapopulation.  In particular, this will get called in two pieces,
// with the seed output collected up for general dispersal.
template <class Individual>
void Patch<Individual>::step_stochastic() {
  deaths();
  add_seeds(births());
}

template <class Individual>
void Patch<Individual>::deaths() {
  for ( species_iterator sp = species.begin();
	sp != species.end(); sp++ )
    sp->deaths();
}

template <class Individual>
std::vector<int> Patch<Individual>::births() {
  std::vector<int> ret(size(), 0);
  std::vector<int>::iterator n = ret.begin();
  for ( species_iterator sp = species.begin();
	sp != species.end(); sp++ )
    *n = sp->births();
  return ret;
}

template <class Individual>
void Patch<Individual>::add_seeds(std::vector<int> seeds) {
  seeds = germination(seeds);
  add_seedlings(seeds);
}

template <class Individual>
void Patch<Individual>::add_seedlings(std::vector<int> seeds) {
  for ( size_t i = 0; i < seeds.size(); i++ )
    species[i].add_seeds(seeds[i]);
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
  for ( size_t i = 0; i < seeds.size(); i++ ) {
    if ( seeds[i] > 0 ) {
      const double p = p_dispersal *
	species[i].germination_probability(environment);
      seeds[i] = p > 0 ? (int)Rf_rbinom(seeds[i], p) : 0.0;
    }
  }
  return seeds;
}

// * ODE interface
template <class Individual>
void Patch<Individual>::derivs(double time,
			       ode::iter_const y, ode::iter dydt) {
  ode_values_set(y);
  ode_rates(dydt);
}

template <class Individual>
size_t Patch<Individual>::ode_size() const {
  return ode::ode_size(species.begin(), species.end());
}

// NOTE: In theory, this is only necessary if no variables have
// changed.  This will often be the case on the first call (because of
// the way that derivs() and ode_set_values works, we take the values
// from the model, set them in the ODE solver, then try to re-set them
// in the model.  Obviously on the first use nothing has changed so we
// should not bother doing anything hard like computing the light
// environment.
template <class Individual>
ode::iter_const Patch<Individual>::ode_values_set(ode::iter_const it) {
  it = ode::ode_values_set(species.begin(), species.end(), it);
  compute_light_environment();
  compute_vars_phys();
  return it;
}

template <class Individual>
ode::iter Patch<Individual>::ode_values(ode::iter it) const {
  return ode::ode_values(species.begin(), species.end(), it);
}

template <class Individual>
ode::iter Patch<Individual>::ode_rates(ode::iter it) const {
  return ode::ode_rates(species.begin(), species.end(), it);
}

// * Private functions
// Sets the strategy for each species
template <class Individual>
void Patch<Individual>::initialise() {
  species.clear();

  // This feels really ugly.
  for ( std::vector<Strategy>::iterator 
	  it = parameters->strategies.begin();
	it != parameters->strategies.end(); it++ ) {
    Species<Individual> s(&(*it)); // (iterator -> object -> pointer)
    species.push_back(s);
  }
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
  for ( species_const_iterator sp = species.begin();
	sp != species.end(); sp++ )
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
double Patch<Individual>::canopy_openness(double height) {
  double tot = 0.0;
  for ( species_const_iterator sp = species.begin();
	sp != species.end(); sp++ )
    tot += sp->leaf_area_above(height);
  // NOTE: patch_area does not appear in the EBT model formulation.
  return exp(-parameters->c_ext * tot / parameters->patch_area);
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

// Given the light environment, "apply" it to every species so that
// physiological variables are updated.
template <class Individual>
void Patch<Individual>::compute_vars_phys() {
  for ( species_iterator sp = species.begin();
	sp != species.end(); sp++ )
    sp->compute_vars_phys(environment);
}

// * R interface

// Actually public functions for interrogating & modifying

template <class Individual>
Species<Individual> Patch<Individual>::r_at(size_t idx) const {
  return species.at(util::check_bounds_r(idx, size()));
}

template <class Individual>
Rcpp::List Patch<Individual>::r_get_species() const {
  Rcpp::List ret;
  for ( species_const_iterator sp = species.begin();
	sp != species.end(); sp++ )
    ret.push_back(Rcpp::wrap(*sp));
  return ret;
}

template <class Individual>
Environment Patch<Individual>::r_environment() const {
  return environment;
}

// In contrast with add_seeds(), we must build the light environment
// in case it has not yet been constructed.
template <class Individual>
void Patch<Individual>::r_add_seeds(std::vector<int> seeds) {
  util::check_length(seeds.size(), size());
  compute_light_environment();
  add_seeds(seeds);
}

template <class Individual>
void Patch<Individual>::r_add_seedlings(std::vector<int> seeds) {
  util::check_length(seeds.size(), size());
  add_seedlings(seeds);
}

template <class Individual>
void Patch<Individual>::r_step() {
  Rcpp::RNGScope scope;
  step();
}

template <class Individual>
void Patch<Individual>::r_step_stochastic() {
  Rcpp::RNGScope scope;
  step_stochastic();
}

// Wrapper functions for testing
template <class Individual>
std::vector<int> Patch<Individual>::r_germination(std::vector<int> seeds) {
  util::check_length(seeds.size(), size());
  return germination(seeds);
}

template <class Individual>
Rcpp::List Patch<Individual>::r_height() const {
  Rcpp::List ret;
  for ( species_const_iterator sp = species.begin();
	sp != species.end(); sp++ )
    ret.push_back(Rcpp::wrap(sp->r_height()));
  return ret;
}

template <class Individual>
void Patch<Individual>::r_set_height(Rcpp::List x) {
  util::check_length(x.size(), size());
  for ( size_t i = 0; i < size(); i++ )
    species[i].r_set_height(x[i]);
}

template <class Individual>
std::vector<int> Patch<Individual>::r_n_individuals() const {
  std::vector<int> n;
  for ( species_const_iterator sp = species.begin();
	sp != species.end(); sp++ )
    n.push_back(sp->r_n_individuals());
  return n;
}

template <class Individual>
void Patch<Individual>::clear() {
  for ( species_iterator sp = species.begin();
	sp != species.end(); sp++ )
    sp->clear();
  environment.clear();
  ode_solver.reset();
}

}

// NOTE: I've not chased up why, but I apparently need to use
// RCPP_EXPOSED_CLASS_NOECL here (because these are templated, not
// real), rather than RCPP_EXPOSED_CLASS.
RCPP_EXPOSED_CLASS_NODECL(model::Patch<model::Plant>)
RCPP_EXPOSED_CLASS_NODECL(model::Patch<model::CohortDiscrete>)

#endif
