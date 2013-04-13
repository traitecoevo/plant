#include "plant_spline.h"
#include "util.h"

namespace model {

PlantSpline::PlantSpline(Strategy s, int n_plants)
  : standalone(true),
    strategy(new Strategy(s)),
    seed(strategy),
    plants_approx(ode_size()) {
  initialise(n_plants);
}

PlantSpline::PlantSpline(Strategy *s, int n_plants)
  : standalone(false),
    strategy(s),
    seed(strategy),
    plants_approx(ode_size()) {
  initialise(n_plants);
}

PlantSpline::PlantSpline(const PlantSpline &other)
  : standalone(other.standalone),
    strategy(standalone ? new Strategy(*other.strategy) : other.strategy),
    seed(other.seed),
    mass_leaf(other.mass_leaf),
    plants(other.plants),
    plants_approx(other.plants_approx) {
}

PlantSpline& PlantSpline::operator=(PlantSpline rhs) {
  swap(*this, rhs);
  return *this;
}

PlantSpline::~PlantSpline() {
  if ( standalone )
    delete strategy;
}

void swap(PlantSpline &a, PlantSpline &b) {
  using std::swap;
  swap(a.standalone,    b.standalone);
  swap(a.strategy,      b.strategy);
  swap(a.seed,          b.seed);
  swap(a.mass_leaf,     b.mass_leaf);
  swap(a.plants,        b.plants);
  swap(a.plants_approx, b.plants_approx);
}

double PlantSpline::max_mass_leaf() const {
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
  if ( m > max_mass_leaf() )
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

// TODO: Get mass_leaf boundaries here; the minimum and maximum
// possible (i.e., size at birth, size at maturity).  This requires
// some poking in Plant and/or Strategy.  For now using hardcoded
// values.  If need be, we can run this out using the ODE step, but
// easier will be to set mass_leaf and look for growth rates heading
// to zero.

// TODO: may want to log transform the spacing of points, or even the
// full spline basis, but not sure.  Do we care more about accuracy at
// the bottom end than the top?
void PlantSpline::initialise(int n_plants) {
  if ( n_plants < 5 )
    ::Rf_error("Need at least 5 plants");

  const double mass_leaf_min = seed.get_mass_leaf();
  // TODO: Really ugly; just guessed and hard coded for now.
  const double mass_leaf_max = 5;

  mass_leaf = util::seq_len(mass_leaf_min, mass_leaf_max, n_plants);

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
