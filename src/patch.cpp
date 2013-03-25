#include "util.h"

#include "functor.h"
#include "patch.h"

namespace model {

Patch::Patch(Parameters p)
  : standalone(true),
    parameters(new Parameters(p)) {
  set_strategies();
}

Patch::Patch(Parameters *p)
  : standalone(false),
    parameters(p) {
  set_strategies();
}

Patch::Patch(const Patch &other)
  : standalone(other.standalone) {
  Rprintf("Copy constructor\n");
  if ( standalone )
    parameters = new Parameters(*other.parameters);
  else
    parameters = other.parameters;
  set_strategies();
}

Patch& Patch::operator=(const Patch &rhs) {
  Rprintf("Assigmnent operator\n");
  // TODO: Violates DRY - must be some way of doing both.
  standalone = rhs.standalone;
  if ( standalone )
    parameters = new Parameters(*rhs.parameters);
  else
    parameters = rhs.parameters;
  set_strategies();

  return *this;
}

Patch::~Patch() {
  if ( standalone )
    delete parameters;
}

size_t Patch::size() const {
  return species.size();
}

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

void Patch::compute_light_environment() {
  // Naive version -- push out to to body of class
  util::Functor<Patch, &Patch::canopy_openness> fun(this);
  light_environment.set_bounds(0, height_max());
  light_environment.set_target(&fun);
  // TODO: should be construct(&fun, 0, height())
  light_environment.construct_spline();
}

spline::Spline Patch::get_light_environment() const {
  return light_environment;
}

void Patch::compute_vars_phys() {
  for ( std::vector<Species>::iterator sp = species.begin();
	sp != species.end(); sp++ )
    sp->compute_vars_phys(&light_environment);
}

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


Rcpp::List Patch::get_plants(int idx) const {
  util::check_bounds(idx, size());
  return species[idx].get_plants();
}

void Patch::add_seed(int idx) {
  species[idx].add_seed();
}

void Patch::r_add_seed(int idx) {
  util::check_bounds(idx, size());
  add_seed(idx);
}

std::vector<double> Patch::r_get_mass_leaf(int idx) const {
  util::check_bounds(idx, size());
  return species[idx].r_get_mass_leaf();
}

void Patch::r_set_mass_leaf(std::vector<double> x, int idx) {
  util::check_bounds(idx, size());
  species[idx].r_set_mass_leaf(x);
}

void Patch::set_strategies() {
  species.clear();

  // This is really ugly.
  for ( std::vector<Strategy>::iterator 
	  it = parameters->strategies.begin();
	it != parameters->strategies.end(); it++ ) {
    Species s(&(*it)); // ugly (iterator -> object -> pointer)
    species.push_back(s);
  }
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

std::vector<double> Patch::r_derivs(std::vector<double> y) {
  // TODO: this check into utils.
  if ( y.size() != ode_size() )
    Rf_error("Incorrect size input (expected %d, recieved %d)",
	     ode_size(), y.size());
  std::vector<double> dydt(y.size());
  const double time = 0.0; // always zero here.
  derivs(time, y.begin(), dydt.begin());
  return dydt;
}

void Patch::r_ode_values_set(std::vector<double> y) {
  // TODO: this check into utils.
  if ( y.size() != ode_size() )
    Rf_error("Incorrect size input (expected %d, recieved %d)",
	     ode_size(), y.size());
  bool changed = false;
  ode_values_set(y.begin(), changed);
}

std::vector<double> Patch::r_ode_values() const {
  std::vector<double> values(ode_size());
  ode_values(values.begin());
  return values;
}

std::vector<double> Patch::r_ode_rates() const {
  std::vector<double> rates(ode_size());
  ode_rates(rates.begin());
  return rates;
}

}
