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

private:
  Strategy *strategy; // do we actually need this?
  std::list<Plant> plants;
};

}

#endif
