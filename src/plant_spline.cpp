#include "plant_spline.h"
#include "util.h"

namespace model {

PlantSpline::PlantSpline(Strategy s, double mass_leaf_max, int n_plants)
  : WithStrategy(s),
    seed(strategy),
    plants_approx(ode_size()) {
  initialise(mass_leaf_max, n_plants);
}

PlantSpline::PlantSpline(Strategy *s, double mass_leaf_max, int n_plants)
  : WithStrategy(s),
    seed(strategy),
    plants_approx(ode_size()) {
  initialise(mass_leaf_max, n_plants);
}

double PlantSpline::mass_leaf_max() const {
  return mass_leaf.back();
}

void PlantSpline::compute_vars_phys(spline::Spline *env) {
  for ( std::vector<Plant>::iterator p = plants.begin();
	p != plants.end(); p++ )
    p->compute_vars_phys(env);
  build_plants_approx();
}

// TODO: This would be a bit nicer if the iterator control was in
// the MultiSpline, but that requires learning about
// ForwardIterator, etc.
ode::iter PlantSpline::ode_rates(double m, ode::iter it) const {
  if ( m > mass_leaf_max() )
    ::Rf_error("Requested plant too large");
  for ( size_t i = 0; i < ode_size(); i++ )
    *it++ = plants_approx.eval(m, i);
  return it;
}

// * R interface
void PlantSpline::r_compute_vars_phys(spline::Spline env) {
  compute_vars_phys(&env);
}

std::vector<double> PlantSpline::r_ode_rates(double m) const {
  std::vector<double> ret(ode_size());
  ode_rates(m, ret.begin());
  return ret;
}

Rcpp::List PlantSpline::r_get_plants() const {
  Rcpp::List ret;
  for ( std::vector<Plant>::const_iterator p = plants.begin();
	p != plants.end(); p++ )
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
// mass_leaf here?  It looks to me that growth in leaf mass only
// asymptotically approaches zero, so there is no value truely large
// enough.  Better to have the really large plants done analytically
// anyway.
void PlantSpline::initialise(double mass_leaf_max, int n_plants) {
  if ( n_plants < 5 )
    ::Rf_error("Need at least 5 plants");

  mass_leaf = util::seq_len(seed.get_mass_leaf(), mass_leaf_max, n_plants);

  Plant p(seed);
  plants.clear(); // defensive, as only used in constructors.
  for ( std::vector<double>::iterator m = mass_leaf.begin();
	m != mass_leaf.end(); m++ ) {
    p.set_mass_leaf(*m);
    plants.push_back(p);
  }

  // Build the approximate plants:
  build_plants_approx();
}

// NOTE: An alternative way would be to clear out the y values only
// (an update_y method, taking an index?).  Could avoid some copying.
void PlantSpline::build_plants_approx() {
  plants_approx.reset();

  std::vector<double> ode_rates_p(plants.begin()->ode_size());

  std::vector<double>::const_iterator m = mass_leaf.begin();
  std::vector<Plant>::iterator p = plants.begin();

  while ( p != plants.end() ) {
    p->ode_rates(ode_rates_p.begin());
    plants_approx.add_point(*m, ode_rates_p);
    ++p;
    ++m;
  }

  plants_approx.init_self();
}

}
