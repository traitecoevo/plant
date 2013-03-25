#include "species.h"

#include "util.h" // is_decreasing

namespace model {

Species::Species() : strategy(NULL) { 
}

Species::Species(Strategy *s) : strategy(s) {
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

void Species::add_seed() {
  Plant p(strategy);
  plants.push_back(p);
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

std::vector<double>::const_iterator
Species::set_values(std::vector<double>::const_iterator it, 
		    bool &changed) {
  for ( std::list<Plant>::iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->set_values(it, changed);
  return it;
}

std::vector<double>::iterator
Species::get_values(std::vector<double>::iterator it) const {
  for ( std::list<Plant>::const_iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->get_values(it);
  return it;
}

std::vector<double>::iterator
Species::get_rates(std::vector<double>::iterator it) const {
  for ( std::list<Plant>::const_iterator p = plants.begin();
	p != plants.end(); p++ )
    it = p->get_rates(it);
  return it;
}



}
