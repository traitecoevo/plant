// -*-c++-*-
#ifndef TREE_SPECIES_
#define TREE_SPECIES_

#include <list>

#include "ode_target.h"
#include "ode_solver.h"
#include "spline.h"
#include "strategy.h"
#include "plant.h"

namespace model {

class Species : public ode::OdeTarget {
public:
  Species();
  Species(Strategy *s);

  // Basic interrogation
  size_t size() const;
  double height_max() const;
  double leaf_area_above(double height) const;

  void compute_vars_phys(spline::Spline *light_environment);

  // ODE interface
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it, bool &changed);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * Births
  int offspring();
  
  // * Deaths
  bool died();

  // Try and get everything.
  Rcpp::List get_plants() const;

  void add_seeds(int n);

  // This is likely to change as more is written, as we'll have to set
  // leaf mass directly anyway.  However, the check within will
  // remain.
  std::vector<double> r_get_mass_leaf() const;
  void r_set_mass_leaf(std::vector<double> x);

private:
  Strategy *strategy;
  Plant seed;
  std::list<Plant> plants;
};

}

#endif
