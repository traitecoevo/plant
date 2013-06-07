// -*-c++-*-
#ifndef TREE_COHORT_TOP_H_
#define TREE_COHORT_TOP_H_

#include "ode_target.h"
#include "strategy.h"
#include "plant.h"

namespace model {

class CohortTop : public Plant {
public:
  CohortTop(Strategy s);
  CohortTop(Strategy *s);

  void compute_vars_phys(const Environment& environment);
  void compute_initial_conditions(const Environment& environment);

  double leaf_area_above(double z) const;

  int offspring();
  bool died();

  // * ODE interface
  size_t ode_size() const;
  ode::iter_const set_ode_values(double time, ode::iter_const it);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * R interface
  double r_growth_rate_gradient(const Environment& environment) const;

private:
  double growth_rate_gradient(const Environment& environment) const;
  double growth_rate_given_height(double height,
				  const Environment& environment);

  double density;
  double density_rate;
  double seeds_survival_weighted;
  double seeds_survival_weighted_rate;

  double time_of_birth;

  static const int ode_dimension = 4;
};

// Painfully copied from Plant, but because we need a non-const
// version.
template <class T, class T2, double (T::*target)(double, T2)>
class FunctorBind2 : public util::DFunctor {
public:
  FunctorBind2(T *obj, T2 arg2) : obj(obj), arg2(arg2) {}
  virtual double operator()(double x) {
    return (obj->*target)(x, arg2);
  }
private:
  T* obj;
  T2 arg2;
};

}

RCPP_EXPOSED_CLASS(model::CohortTop)

#endif
