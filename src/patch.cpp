#include "util.h"

#include "functor.h"
#include "patch.h"

namespace model {

Patch::Patch(Parameters p)
  : standalone(true),
    parameters(new Parameters(p)),
    ode_solver(this) {
  initialise();
}

Patch::Patch(Parameters *p)
  : standalone(false),
    parameters(p),
    ode_solver(this) {
  initialise();
}

Patch::Patch(const Patch &other)
  : standalone(other.standalone), 
    ode_solver(this) {
  Rprintf("Copy constructor\n");
  if ( standalone )
    parameters = new Parameters(*other.parameters);
  else
    parameters = other.parameters;
  initialise();
}

Patch& Patch::operator=(Patch rhs) {
  swap(*this, rhs);
  return *this;
}

Patch::~Patch() {
  if ( standalone )
    delete parameters;
}

// void Patch::step() {}

// TODO: This is a hack for now, as the state setter is using the R
// function.  Instead set_state should at least offer to take an
// interator.
// 
// TODO: This should only move in state that needs changing, or update
// things if the dimension has changed.  For now that, requires
// knowledge about how I'm going to do things I've not decided yet.
void Patch::step_deterministic() {
  double time = 0.0; // TODO: Move up?  Ignore?
  // TODO: Oh dear -- relying on an R-only function here.
  ode_solver.set_state(r_ode_values(), time);
  ode_solver.step();
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
    Species s(&(*it)); // ugly (iterator -> object -> pointer)
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

std::vector<double> Patch::r_get_mass_leaf(int idx) const {
  util::check_bounds(idx, size());
  return species[idx].r_get_mass_leaf();
}

void Patch::r_set_mass_leaf(std::vector<double> x, int idx) {
  util::check_bounds(idx, size());
  species[idx].r_set_mass_leaf(x);
}

}
