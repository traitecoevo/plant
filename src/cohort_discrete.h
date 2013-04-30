// -*-c++-*-
#ifndef TREE_COHORT_DISCRETE_H_
#define TREE_COHORT_DISCRETE_H_

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

  // * R interface
  int r_n_individuals() const;
  void r_set_n_individuals(int n);

private:
  int n_individuals;
};

}

RCPP_EXPOSED_CLASS(model::CohortDiscrete)

#endif
