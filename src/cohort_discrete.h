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
  CohortDiscrete(Strategy  s, int n_individuals_);
  CohortDiscrete(Strategy *s, int n_individuals_);

  double leaf_area_above(double z) const;
  int offspring();
  bool died();

  // * R interface
  int get_n_individuals() const;
  void set_n_individuals(int n);

  size_t state_size() const;
  state::iterator       get_state(state::iterator       it) const;
  state::const_iterator set_state(state::const_iterator it);

private:
  int n_individuals;
};

}

RCPP_EXPOSED_CLASS(model::CohortDiscrete)

#endif
