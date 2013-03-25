#include "species.h"

namespace model {

Species::Species() : strategy(NULL) { 
}

Species::Species(Strategy *s) : strategy(s) {
}

size_t Species::size() const {
  return plants.size();
}

double Species::height_max() const {
  return plants.begin()->get_height();
}

double Species::leaf_area_above(double height) const {
  double tot = 0.0;
  for ( std::list<Plant>::const_iterator it = plants.begin();
	it != plants.end(); it++ )
    tot += it->leaf_area_above(height);
  return tot;
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

}
