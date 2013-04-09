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

// This class is only used by Patch, and will possibly be where some
// major changes occur once we move to allowing Cohorts.

// It's possible that we should just friend this with Patch.

class Species : public ode::OdeTarget {
public:
  Species();
  Species(Strategy *s);

  // * Births and deaths
  int births();
  int deaths();

  // * Lower level functions, used by Patch
  size_t size() const;
  double height_max() const;
  double leaf_area_above(double height) const;
  void compute_vars_phys(spline::Spline *light_environment);
  void add_seeds(int n);
  void clear();

  // * ODE interface
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it, bool &changed);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * R interface
  // Even though there is not an R interface to this class, these
  // functions are used only by other R interface functions.
  std::vector<double> r_get_mass_leaf() const;
  void r_set_mass_leaf(std::vector<double> x);
  Rcpp::List r_get_plants() const;

private:
  Strategy *strategy;
  Plant seed;
  std::list<Plant> plants;
};

}

#endif
