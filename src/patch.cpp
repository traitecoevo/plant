#include "util.h"

#include "functor.h"
#include "patch.h"

namespace model {

Patch::Patch(Parameters p)
  : standalone(true),
    parameters(new Parameters(p)),
    age(0.0),
    ode_solver(this) {
  initialise();
}

Patch::Patch(Parameters *p)
  : standalone(false),
    parameters(p),
    age(0.0),
    ode_solver(this) {
  initialise();
}

Patch::Patch(const Patch &other)
  : standalone(other.standalone),
    parameters(standalone ?
	       new Parameters(*other.parameters) : other.parameters),
    age(other.age),
    light_environment(other.light_environment),
    species(other.species),
    ode_solver(this) {
}

Patch& Patch::operator=(Patch rhs) {
  swap(*this, rhs);
  return *this;
}

Patch::~Patch() {
  if ( standalone )
    delete parameters;
}

void Patch::step() {
  step_deterministic();
  step_stochastic();
}

void Patch::step_deterministic() {
  std::vector<double> y(ode_size());
  ode_values(y.begin());
  ode_solver.set_state(y, age);
  ode_solver.step();
  age = ode_solver.get_time();
}

// TODO: Should this be a method within species, perhaps?  If so then,
// we get to do a std::foreach.  Probably not, because that means that
// the Species class contains code depending on how it is used (wrong
// level).
void Patch::step_stochastic() {
  for ( std::vector<Species>::iterator sp = species.begin();
	sp != species.end(); sp++ ) {
    sp->deaths();
    sp->add_seeds(sp->births());
  }
}

// * ODE interface
void Patch::derivs(double time,
		   std::vector<double>::const_iterator y,
		   std::vector<double>::iterator dydt) {
  bool changed = false;
  ode_values_set(y, changed);
  // Next two will be optional (if (changed))
  compute_light_environment();
  compute_vars_phys();

  ode_rates(dydt);
}

size_t Patch::ode_size() const {
  size_t ret = 0;
  for ( std::vector<Species>::const_iterator sp = species.begin();
	sp != species.end(); sp++ )
    ret += sp->ode_size();
  return ret;
}

ode::iter_const Patch::ode_values_set(ode::iter_const it, bool &changed) {
  for ( std::vector<Species>::iterator sp = species.begin();
	sp != species.end(); sp++ )
    it = sp->ode_values_set(it, changed);
  return it;
}

ode::iter Patch::ode_values(ode::iter it) const {
  for ( std::vector<Species>::const_iterator sp = species.begin(); 
	sp != species.end(); sp++ )
    it = sp->ode_values(it);
  return it;
}

ode::iter Patch::ode_rates(ode::iter it) const {
  for ( std::vector<Species>::const_iterator sp = species.begin(); 
	sp != species.end(); sp++ )
    it = sp->ode_rates(it);
  return it;
}

// * Private functions
void Patch::swap(Patch &a, Patch &b) {
  using std::swap;
  swap(a.standalone,        b.standalone);
  swap(a.parameters,        b.parameters);
  swap(a.age,               b.age);
  swap(a.light_environment, b.light_environment);
  swap(a.species,           b.species);
  swap(a.ode_solver,        b.ode_solver);
}

// Sets the strategy for each species
void Patch::initialise() {
  species.clear();

  // This feels really ugly.
  for ( std::vector<Strategy>::iterator 
	  it = parameters->strategies.begin();
	it != parameters->strategies.end(); it++ ) {
    Species s(&(*it)); // (iterator -> object -> pointer)
    species.push_back(s);
  }
}

// Number of species
size_t Patch::size() const {
  return species.size();
}

// Maxiumum height for any species in the Patch.  Empty patches (no
// species or no individuals) have height 0.
double Patch::height_max() const {
  double ret = 0.0;
  for ( std::vector<Species>::const_iterator sp = species.begin();
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
double Patch::canopy_openness(double height) {
  double tot = 0.0;
  for ( std::vector<Species>::const_iterator sp = species.begin();
	sp != species.end(); sp++ )
    tot += sp->leaf_area_above(height);
  // NOTE: patch_area does not appear in the EBT model formulation.
  return exp(-parameters->c_ext * tot / parameters->patch_area);
}

// Create the spline the characterises the light environment.
void Patch::compute_light_environment() {
  // Naive version -- push out to to body of class
  util::Functor<Patch, &Patch::canopy_openness> fun(this);
  light_environment.set_bounds(0, height_max());
  light_environment.set_target(&fun);
  // TODO: should be construct(&fun, 0, height())
  light_environment.construct_spline();
}

// Given the light environment, "apply" it to every species so that
// physiological variables are updated.
void Patch::compute_vars_phys() {
  for ( std::vector<Species>::iterator sp = species.begin();
	sp != species.end(); sp++ )
    sp->compute_vars_phys(&light_environment);
}

// * R interface

// Actually public functions for interrogating & modifying
Rcpp::List Patch::r_get_plants(int idx) const {
  util::check_bounds(idx, size());
  return species[idx].r_get_plants();
}

spline::Spline Patch::r_light_environment() const {
  return light_environment;
}

void Patch::r_add_seed(int idx) {
  util::check_bounds(idx, size());
  species[idx].add_seeds(1);
}

void Patch::r_step() {
  Rcpp::RNGScope scope;
  step();
}

void Patch::r_step_stochastic() {
  Rcpp::RNGScope scope;
  step_stochastic();
}

// Wrapper functions for testing
size_t Patch::r_size() const {
  return size();
}

double Patch::r_height_max() const {
  return height_max();
}

double Patch::r_canopy_openness(double height) {
  return canopy_openness(height);
}

void Patch::r_compute_light_environment() {
  compute_light_environment();
}

void Patch::r_compute_vars_phys() {
  compute_vars_phys();
}

double Patch::r_age() const {
  return age;
}

std::vector<double> Patch::r_get_mass_leaf(int idx) const {
  util::check_bounds(idx, size());
  return species[idx].r_get_mass_leaf();
}

void Patch::r_set_mass_leaf(std::vector<double> x, int idx) {
  util::check_bounds(idx, size());
  species[idx].r_set_mass_leaf(x);
}

}
