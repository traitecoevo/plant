// -*-c++-*-
#ifndef TREE_COHORT_
#define TREE_COHORT_

// #include <Rcpp.h>

#include "ode_target.h"
#include "strategy.h"
#include "plant.h"

namespace model {

class Cohort : public ode::OdeTarget {
public:
  Cohort(Strategy s);
  Cohort(Strategy *s);

  // Copy constructor, assigment and destructor (rule of three)
  Cohort(const Cohort &other);
  Cohort& operator=(Cohort rhs);
  ~Cohort();
  // and a half.
  friend void swap(Cohort &a, Cohort &b);

  // TODO: which height should this return?  Top I think?
  double get_height() const;

  // TODO: which height should this set?
  // [eqn 1-8] Update size variables to a new leaf mass.
  void set_mass_leaf(double mass_leaf_mean, double mass_leaf_top);

  // * ODE interface
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it, bool &changed);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

private:
  // As with other similar components; involved in memory management.
  bool standalone;
  Strategy *strategy;

  // Two plants; one representing the mean individual and the other
  // representing the top.
  Plant plant_mean, plant_top;

  static const int ode_dimension = 7; // I think.
};

}

#endif
