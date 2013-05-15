// -*-c++-*-
#ifndef TREE_COHORT_TOP_H_
#define TREE_COHORT_TOP_H_

#include "ode_target.h"
#include "strategy.h"
#include "plant.h"
#include "disturbance.h"

namespace model {

class CohortTop : public Plant {
public:
  CohortTop(Strategy s);
  CohortTop(Strategy *s);

  // TODO: See design.md (search: compute_vars_phys_surv)
  void compute_vars_phys_surv(spline::Spline *env, double survival_patch);
  void compute_initial_conditions(spline::Spline *env,
				  double seed_input);

  double leaf_area_above(double z) const;
  int offspring();
  bool died();

  // * ODE interface
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * R interface
  void r_compute_vars_phys_surv(spline::Spline env,
				double survival_patch);
  void r_compute_initial_conditions(spline::Spline env,
				    double seed_input);
  double r_growth_rate_gradient(spline::Spline env) const;
  double r_growth_rate_given_mass(double mass_leaf,
				  spline::Spline env);

private:
  double growth_rate_gradient(spline::Spline *env) const;
  double growth_rate_given_mass(double mass_leaf, spline::Spline *env);

  double density;
  double density_rate;
  double seeds_survival_weighted;
  double seeds_survival_weighted_rate;

  double age_at_birth;

  static const int ode_dimension = 4;
};

// Painfully copied from Plant, but because we need a non-const
// version.
//
// TODO: Might be a solution here.
// http://stackoverflow.com/questions/14926482/const-and-non-const-template-specialization
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
