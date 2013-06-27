#include "plant_spline.h"
#include "util.h"

namespace model {

PlantSpline::PlantSpline(Strategy s, double height_max, int n_plants)
  : strategy(s),
    seed(strategy.get()),
    plants_approx(ode_size()) {
  initialise(height_max, n_plants);
}

PlantSpline::PlantSpline(Strategy *s, double height_max, int n_plants)
  : strategy(s),
    seed(strategy.get()),
    plants_approx(ode_size()) {
  initialise(height_max, n_plants);
}

double PlantSpline::height_max() const {
  return height.back();
}

void PlantSpline::compute_vars_phys(const Environment& environment) {
  for (std::vector<Plant>::iterator p = plants.begin();
       p != plants.end(); ++p)
    p->compute_vars_phys(environment);
  build_plants_approx();
}

ode::iter PlantSpline::ode_rates(double height, ode::iter it) const {
  if ( height > height_max() )
    ::Rf_error("Requested plant too large");
  for ( size_t i = 0; i < ode_size(); i++ )
    *it++ = plants_approx.eval(height, i);
  return it;
}

Rcpp::NumericVector PlantSpline::r_get_vars_phys(double m) const {
  Rcpp::NumericVector ret = seed.r_get_vars_phys();
  for ( size_t i = 0; i < ret.size(); i++ )
    ret[i] = NA_REAL;
  return ret;
}

std::vector<double> PlantSpline::r_ode_rates(double height) const {
  std::vector<double> ret(ode_size());
  ode_rates(height, ret.begin());
  return ret;
}

Rcpp::List PlantSpline::r_get_plants() const {
  Rcpp::List ret;
  for (std::vector<Plant>::const_iterator p = plants.begin();
       p != plants.end(); ++p)
    ret.push_back(*p);
  return ret;
}

spline::MultiSpline PlantSpline::r_get_plants_approx() const {
  return plants_approx;
}

// * Private methods
size_t PlantSpline::ode_size() const {
  return seed.ode_size();
}

// NOTE: Perhaps more gracefully deal with the upper boundary of
// height here?  It looks to me that growth in height only
// asymptotically approaches zero, so there is no value truely large
// enough.  Better to have the really large plants done analytically
// anyway.
void PlantSpline::initialise(double height_max, int n_plants) {
  if ( n_plants < 5 )
    ::Rf_error("Need at least 5 plants");

  height = util::seq_len(seed.height(), height_max, n_plants);

  Plant p(seed);
  plants.clear(); // defensive, as only used in constructors.
  for (std::vector<double>::iterator h = height.begin();
       h != height.end(); ++h) {
    p.set_height(*h);
    plants.push_back(p);
  }

  // Build the approximate plants:
  build_plants_approx();
}

// NOTE: An alternative way would be to clear out the y values only
// (an update_y method, taking an index?).  Could avoid some copying.
void PlantSpline::build_plants_approx() {
  plants_approx.clear();

  std::vector<double> ode_rates_p(plants.begin()->ode_size());

  std::vector<double>::const_iterator h = height.begin();
  std::vector<Plant>::iterator p = plants.begin();

  while ( p != plants.end() ) {
    p->ode_rates(ode_rates_p.begin());
    plants_approx.add_point(*h, ode_rates_p);
    ++p;
    ++h;
  }

  plants_approx.init_self();
}

}
