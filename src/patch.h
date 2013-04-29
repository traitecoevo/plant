// -*-c++-*-
#ifndef TREE_PATCH_
#define TREE_PATCH_

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
  virtual size_t r_size() const = 0;
  virtual double r_height_max() const = 0;
  virtual double r_canopy_openness(double) = 0;
  virtual void r_compute_light_environment() = 0;
  virtual spline::Spline r_light_environment() const = 0;
  virtual void r_compute_vars_phys() = 0;
  virtual double r_age() const = 0;
  virtual std::vector<int> r_germination(std::vector<int> seeds) = 0;
  virtual Rcpp::List r_get_species() const = 0;
  virtual void r_add_seeds(std::vector<int> seeds) = 0;
  virtual void r_add_seedlings(std::vector<int> seeds) = 0;
  virtual std::vector<double> r_get_mass_leaf(size_t idx) const = 0;
  virtual void r_set_mass_leaf(std::vector<double> x, size_t idx) = 0;
  virtual std::vector<int> r_n_individuals() const = 0;
  virtual void r_clear() = 0;
  virtual void r_step() = 0;
  virtual void step_deterministic() = 0;
  virtual void r_step_stochastic() = 0;
};

template <class Individual>
class Patch : public PatchBase {
public:
  Patch(Parameters p);
  Patch(Parameters *p);

  // Advance the system through one complete time step.
  void step();

  // Advance the system through one time step deterministically
  // (plant growth, physiological accounting)
  void step_deterministic();
  // Advance the system through the stochastic life cycle stages
  // (producing seeds and dying).
  void step_stochastic();
  // TODO: Not sure (until Metapopulation developed) if these should
  // be public or private.
  std::vector<int> births();
  void deaths();
  std::vector<int> germination(std::vector<int> seeds);
  void add_seeds(std::vector<int> seeds);
  void add_seedlings(std::vector<int> seeds);

  // * ODE interface.
  void derivs(double time, ode::iter_const y, ode::iter dydt);
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * R interface.

  // Actually public functions for interrogating & modifying
  Species<Individual> r_at(size_t idx) const;
  Rcpp::List r_get_species() const;
  spline::Spline r_light_environment() const;
  void r_add_seeds(std::vector<int> seeds);
  void r_add_seedlings(std::vector<int> seeds);
  void r_step();
  void r_step_stochastic(); // step_stochastic, plus RNG control

  // Wrapper functions for testing.  Some of these expose otherwise
  // private methods, and do nothing more than pass through.
  size_t r_size() const { return size(); }
  double r_height_max() const { return height_max(); }
  double r_canopy_openness(double height) {return canopy_openness(height);}
  void r_compute_light_environment() {compute_light_environment();}
  void r_compute_vars_phys() {compute_vars_phys();}
  double r_age() const {return age;}

  // These include size or bounds checking
  std::vector<int> r_germination(std::vector<int> seeds);
  // TODO: This is likely to change as more is written.
  std::vector<double> r_get_mass_leaf(size_t idx) const;
  void r_set_mass_leaf(std::vector<double> x, size_t idx);
  std::vector<int> r_n_individuals() const;
  
  // TODO: Should be clear()?
  void r_clear();
  
private:
  void initialise();

  // Number of species
  size_t size() const;

  // Maximum height for any species in the Patch
  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  double canopy_openness(double height);

  void compute_light_environment();
  void compute_vars_phys();

  Parameters::ptr parameters;
  double age;

  spline::AdaptiveSpline light_environment;

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
    age(0.0),
    ode_solver(this) {
  initialise();
}

template <class Individual>
Patch<Individual>::Patch(Parameters *p)
  : parameters(p),
    age(0.0),
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
  ode_solver.set_state(y, age);
  ode_solver.step();
  age = ode_solver.get_time();
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
      const double p =
	p_dispersal * species[i].germination_probability(&light_environment);
      if ( p > 0 )
	seeds[i] = (int)Rf_rbinom(seeds[i], p);
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
  size_t ret = 0;
  for ( species_const_iterator sp = species.begin();
	sp != species.end(); sp++ )
    ret += sp->ode_size();
  return ret;
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
  for ( species_iterator sp = species.begin();
	sp != species.end(); sp++ )
    it = sp->ode_values_set(it);

  compute_light_environment();
  compute_vars_phys();

  return it;
}

template <class Individual>
ode::iter Patch<Individual>::ode_values(ode::iter it) const {
  for ( species_const_iterator sp = species.begin(); 
	sp != species.end(); sp++ )
    it = sp->ode_values(it);
  return it;
}

template <class Individual>
ode::iter Patch<Individual>::ode_rates(ode::iter it) const {
  for ( species_const_iterator sp = species.begin(); 
	sp != species.end(); sp++ )
    it = sp->ode_rates(it);
  return it;
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

// Maxiumum height for any species in the Patch.  Empty patches (no
// species or no individuals) have height 0.
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
template <class Individual>
void Patch<Individual>::compute_light_environment() {
  // Naive version -- push out to to body of class
  util::Functor<Patch, &Patch<Individual>::canopy_openness> fun(this);
  light_environment.set_bounds(0, height_max());
  light_environment.set_target(&fun);
  // TODO: should be construct(&fun, 0, height())
  light_environment.construct_spline();
}

// Given the light environment, "apply" it to every species so that
// physiological variables are updated.
template <class Individual>
void Patch<Individual>::compute_vars_phys() {
  for ( species_iterator sp = species.begin();
	sp != species.end(); sp++ )
    sp->compute_vars_phys(&light_environment);
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
spline::Spline Patch<Individual>::r_light_environment() const {
  return light_environment;
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
std::vector<double> Patch<Individual>::r_get_mass_leaf(size_t idx) const {
  util::check_bounds(idx, size());
  return species[idx].r_get_mass_leaf();
}

template <class Individual>
void Patch<Individual>::r_set_mass_leaf(std::vector<double> x, size_t idx) {
  util::check_bounds(idx, size());
  species[idx].r_set_mass_leaf(x);
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
void Patch<Individual>::r_clear() {
  age = 0.0;
  for ( species_iterator sp = species.begin();
	sp != species.end(); sp++ )
    sp->clear();
  ode_solver.reset();
}

}

// NOTE: I've not chased up why, but I apparently need to use
// RCPP_EXPOSED_CLASS_NOECL here (because these are templated, not
// real), rather than RCPP_EXPOSED_CLASS.
RCPP_EXPOSED_CLASS_NODECL(model::Patch<model::Plant>)
RCPP_EXPOSED_CLASS_NODECL(model::Patch<model::CohortDiscrete>)

#endif
