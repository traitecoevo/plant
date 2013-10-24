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

  // NOTE: I'd rather just make these private to disable, but that
  // conflicts with Rcpp's inheritance model.
  int offspring();
  bool died();

  // * ODE interface
  size_t ode_size() const;
  ode::iterator_const set_ode_values(double time, ode::iterator_const it);
  ode::iterator       ode_values(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it)  const;

  // * R interface
  double r_growth_rate_gradient(const Environment& environment);
  double r_growth_rate_given_height(double height_,
				    const Environment& environment);

private:
  double growth_rate_gradient(const Environment& environment,
			      integration::intervals_type intervals) const;
  double growth_rate_given_height(double height_,
				  const Environment& environment);

  double log_density;
  double log_density_rate;
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
  FunctorBind2(T *obj_, T2 arg2_) : obj(obj_), arg2(arg2_) {}
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
