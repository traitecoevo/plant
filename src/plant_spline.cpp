#include "plant_spline.h"
#include "util.h"

namespace model {

PlantSpline::PlantSpline(Strategy s, int n_plants)
  : standalone(true),
    strategy(new Strategy(s)) {
  build(n_plants);
}

PlantSpline::PlantSpline(Strategy *s, int n_plants)
  : standalone(false),
    strategy(s) {
  build(n_plants);
}

PlantSpline::PlantSpline(const PlantSpline &other)
  : standalone(other.standalone),
    strategy(standalone ? new Strategy(*other.strategy) : other.strategy),
    mass_leaf_log(other.mass_leaf_log),
    plants(other.plants),
    ode_values(other.ode_values),
    ode_values_approx(other.ode_values_approx) {
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
  swap(a.standalone,        b.standalone);
  swap(a.strategy,          b.strategy);
  swap(a.mass_leaf_log,     b.mass_leaf_log);
  swap(a.plants,            b.plants);
  swap(a.ode_values,        b.ode_values);
  swap(a.ode_values_approx, b.ode_values_approx);
}

// This is hard and ugly because we want to fill the matrix by row and
// access it by column.  Could do in three passes, but would have a
// bad effect on the amount of copying.
void PlantSpline::compute_vars_phys(spline::Spline *env) {
  const size_t n_var = ode_size();
  std::vector<double> ode_values_p(n_var);

  size_t i = 0;
  for ( std::vector<Plant>::iterator p = plants.begin();
	p != plants.end(); p++, i++ ) {
    // Actually compute the physiological variables.
    p->compute_vars_phys(env);

    // Copy these into a little local vector.
    p->ode_values(ode_values_p.begin());
    
    // And copy these into our big 2d vector:
    for ( int j = 0; j < n_var; j++ )
      ode_values[j][i] = ode_values_p[j];
  }

  // Then build splines:
  for ( size_t j = 0; j < n_var; j++ )
    ode_values_approx[j].init(mass_leaf_log, ode_values[j]);
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

Rcpp::NumericMatrix PlantSpline::r_get_ode_values() const {
  const size_t n_var = ode_size();
  Rcpp::NumericMatrix ret(plants.size(), n_var);
  for ( size_t i = 0; i < n_var; i++ )
    for ( size_t j = 0; j < plants.size(); j++ )
      ret(j,i) = ode_values[i][j];
  return ret;
}

Rcpp::List PlantSpline::r_get_ode_values_approx() const {
  Rcpp::List ret;
  for ( std::vector<spline::Spline>::const_iterator 
	  it = ode_values_approx.begin();
	it != ode_values_approx.end(); it++ )
    ret.push_back(*it);
  return ret;
}

// * Private methods
size_t PlantSpline::ode_size() const {
  if ( plants.size() > 0 ) {
    return plants.begin()->ode_size();
  } else {
    Plant tmp(strategy);
    return tmp.ode_size();
  }
}

void PlantSpline::build(int n_plants) {
  // TODO: Get mass_leaf boundaries here; the minimum and maximum
  // possible (i.e., size at birth, size at maturity).  This requires
  // some poking in Plant and/or Strategy.  For now using hardcoded
  // values.  If need be, we can run this out using the ODE step, but
  // easier will be to set mass_leaf and look for growth rates heading
  // to zero.
  const double mass_leaf_min = 1e-5, mass_leaf_max = 5;
  mass_leaf_log = util::seq_len(log(mass_leaf_min),
				log(mass_leaf_max),
				n_plants);

  plants.clear();
  Plant p(strategy);
  for ( std::vector<double>::iterator it = mass_leaf_log.begin();
	it != mass_leaf_log.end(); it++ ) {
    p.set_mass_leaf(exp(*it));
    plants.push_back(p);
  }

  std::vector<double> tmp(n_plants);
  for ( size_t i = 0; i < ode_size(); i++ )
    ode_values.push_back(tmp);
  ode_values_approx.resize(ode_size());
}

}
