// -*-c++-*-
#ifndef TREE_SPECIES_
#define TREE_SPECIES_

#include <list>

#include "strategy.h"
#include "plant.h"

namespace model {

class Species {
public:
  Species();
  Species(Strategy *s);

  // Basic interrogation
  size_t size() const;
  double height_max() const;
  double leaf_area_above(double height) const;

  // Try and get everything.
  Rcpp::List get_plants() const;

  void add_seed();

  // This is likely to change as more is written, as we'll have to set
  // leaf mass directly anyway.  However, the check within will
  // remain.
  std::vector<double> r_get_mass_leaf() const;
  void r_set_mass_leaf(std::vector<double> x);

private:
  Strategy *strategy;
  std::list<Plant> plants;
};

}

#endif
