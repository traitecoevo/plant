// -*-c++-*-
#ifndef TREE_COHORT_DISCRETE_
#define TREE_COHORT_DISCRETE_

#include "ode_target.h"
#include "strategy.h"
#include "plant.h"

namespace model {

class CohortDiscrete : public ode::OdeTarget {
public:
  CohortDiscrete(Strategy s,  int n_individuals);
  CohortDiscrete(Strategy *s, int n_individuals);

  // Copy constructor, assigment and destructor (rule of three)
  CohortDiscrete(const CohortDiscrete &other);
  CohortDiscrete& operator=(CohortDiscrete rhs);
  ~CohortDiscrete();
  // and a half.
  friend void swap(CohortDiscrete &a, CohortDiscrete &b);

  double get_height() const;
  void set_mass_leaf(double mass_leaf_);
  double leaf_area_above(double z) const;
  void compute_vars_phys(spline::Spline *env);

  int offspring();
  bool died();
  double germination_probability(spline::Spline *env);

  // // * ODE interface
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it, bool &changed);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // // * R interface
  // double r_get_mass_leaf() const;
  // Rcpp::NumericVector r_get_vars_size() const;
  // Rcpp::NumericVector r_get_vars_phys() const;
  // Rcpp::List r_get_parameters() const;
  // void r_compute_vars_phys(spline::Spline env);
  // double r_compute_assimilation(spline::Spline env) const;
  // double r_compute_assimilation_x(double x, spline::Spline env) const;
  // double r_germination_probability(spline::Spline env);
  // std::string r_name() const;

private:
  bool standalone;
  Strategy *strategy;
  Plant plant;
  int n_individuals;
};

}

#endif
