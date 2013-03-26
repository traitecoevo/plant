#include "species.h"

#include "util.h" // is_decreasing

namespace model {

// TODO: I'm a bit wary of the impact of seed(NULL) here, especially
// if a copy constructor is triggered, but also just in general.  It
// would be nice if we could just skip this contructor entirely, but
// it is apparently necessary for something in Patch, I think.
Species::Species() : 
  strategy(NULL),
  seed(strategy) {
}

Species::Species(Strategy *s) : 
  strategy(s),
  seed(strategy) {
}

size_t Species::size() const {
  return plants.size();
}

double Species::height_max() const {
  if ( size() == 0 )
    return 0.0;
  else
    return plants.begin()->get_height();
}

double Species::leaf_area_above(double height) const {
  double tot = 0.0;
  for ( std::list<Plant>::const_iterator it = plants.begin();
	it != plants.end(); it++ )
    tot += it->leaf_area_above(height);
  return tot;
}

void Species::compute_vars_phys(spline::Spline *light_environment) {
  for ( std::list<Plant>::iterator it = plants.begin();
	it != plants.end(); it++ )
    it->compute_vars_phys(light_environment);
}

int Species::offspring() {
  int born = 0;
  for ( std::list<Plant>::iterator it = plants.begin();
	it != plants.end(); it++ )
    born += it->offspring();
  return born;
}

bool Species::died() {
  bool any_died = false;
  std::list<Plant>::iterator it = plants.begin();
  while ( it != plants.end() ) {
    if ( it->died() ) {
      any_died = true;
      it = plants.erase(it); // increments iterator
    } else {
      it++;
    }
  }

  return any_died;
}

void Species::add_seeds(int n) {
  for ( ; n > 0; n-- )
    plants.push_back(seed);
}

Rcpp::List Species::get_plants() const {
  Rcpp::List ret;
  for ( std::list<Plant>::const_iterator it = plants.begin();
	it != plants.end(); it++ )
    ret.push_back(Rcpp::wrap(*it));
  return ret;
}

std::vector<double> Species::r_get_mass_leaf() const { 
  std::vector<double> ret;
  std::list<Plant>::const_iterator p = plants.begin(); 
  while ( p != plants.end() )
    ret.push_back((p++)->get_mass_leaf());
  return ret;
}

// NOTE: Roll back on error is not possible here at present.
void Species::r_set_mass_leaf(std::vector<double> x) {
  if ( x.size() != size() )
    Rf_error("Unexpected size of mass_leaf: expected %d, recieved %d",
	     size(), x.size());
  if ( !util::is_decreasing(x.begin(), x.end()) )
    Rf_error("mass_leaf must be decreasing (ties allowed)");
  std::vector<double>::iterator it = x.begin();
  std::list<Plant>::iterator p = plants.begin();
  while ( p != plants.end() ) {
    p->set_mass_leaf(*it++);
    p++;
  }
}

size_t Species::ode_size() const {
  return size() * Plant::ode_dimension;
}

ode::iter_const Species::ode_values_set(ode::iter_const it, bool &changed) {
  for ( std::list<Plant>::iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->ode_values_set(it, changed);
  return it;
}

ode::iter Species::ode_values(ode::iter it) const {
  for ( std::list<Plant>::const_iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->ode_values(it);
  return it;
}

ode::iter Species::ode_rates(ode::iter it) const {
  for ( std::list<Plant>::const_iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->ode_rates(it);
  return it;
}



}
