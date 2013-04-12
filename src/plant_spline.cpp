#include "plant_spline.h"
#include "util.h"

namespace model {

PlantSpline::PlantSpline(Strategy s, int n_plants)
  : standalone(true),
    strategy(new Strategy(s)),
    seed(strategy),
    plants_approx(seed.ode_size()) {
  initialise(n_plants);
}

PlantSpline::PlantSpline(Strategy *s, int n_plants)
  : standalone(false),
    strategy(s),
    seed(strategy),
    plants_approx(seed.ode_size()) {
  initialise(n_plants);
}

PlantSpline::PlantSpline(const PlantSpline &other)
  : standalone(other.standalone),
    strategy(standalone ? new Strategy(*other.strategy) : other.strategy),
    seed(other.seed),
    mass_leaf_log(other.mass_leaf_log),
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
  swap(a.standalone,     b.standalone);
  swap(a.strategy,       b.strategy);
  swap(a.seed,           b.seed);
  swap(a.mass_leaf_log,  b.mass_leaf_log);
  swap(a.plants,         b.plants);
  swap(a.plants_approx,  b.plants_approx);
}

void PlantSpline::compute_vars_phys(spline::Spline *env) {
  for ( std::vector<Plant>::iterator p = plants.begin();
	p != plants.end(); p++ )
    p->compute_vars_phys(env);
  build_plants_approx();
}

// * R interface
void PlantSpline::r_compute_vars_phys(spline::Spline env) {
  compute_vars_phys(&env);
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

// TODO: Get mass_leaf boundaries here; the minimum and maximum
// possible (i.e., size at birth, size at maturity).  This requires
// some poking in Plant and/or Strategy.  For now using hardcoded
// values.  If need be, we can run this out using the ODE step, but
// easier will be to set mass_leaf and look for growth rates heading
// to zero.

// TODO: not sure about the handling of log/exp here; is it wanted /
// desirable?  Why do we care more about accuracy at the bottom end
// than the top?
void PlantSpline::initialise(int n_plants) {
  if ( n_plants < 5 )
    ::Rf_error("Need at least 5 plants");

  // TODO: Ugly at the moment (using r_ method)
  const double mass_leaf_min = seed.r_get_mass_leaf();
  // TODO: Really ugly; just guessed and hard coded for now.
  const double mass_leaf_max = 5;

  mass_leaf_log = util::seq_len(log(mass_leaf_min),
				log(mass_leaf_max),
				n_plants);

  Plant p(seed);
  plants.clear(); // defensive, as only used in constructors.
  for ( std::vector<double>::iterator lm = mass_leaf_log.begin();
	lm != mass_leaf_log.end(); lm++ ) {
    p.set_mass_leaf(exp(*lm));
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

  std::vector<double>::const_iterator m = mass_leaf_log.begin();
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
