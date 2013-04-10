// -*-c++-*-
#ifndef TREE_COHORT_DISCRETE_
#define TREE_COHORT_DISCRETE_

#include "ode_target.h"
#include "strategy.h"
#include "plant.h"

namespace model {

class CohortDiscrete : public Plant {
public:
  CohortDiscrete(Strategy  s);
  CohortDiscrete(Strategy *s);
  CohortDiscrete(Strategy  s, int n_individuals);
  CohortDiscrete(Strategy *s, int n_individuals);

  double leaf_area_above(double z) const;
  int offspring();
  bool died();

private:
  int n_individuals;
};

}

#endif
