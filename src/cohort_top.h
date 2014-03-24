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
  double leaf_area() const;
  double fecundity() const;

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
  size_t state_size() const;
  state::iterator       get_state(state::iterator       it) const;
  state::const_iterator set_state(state::const_iterator it);

private:
  double growth_rate_gradient(const Environment& environment,
			      integration::intervals_type intervals) const;
  double growth_rate_given_height(double height_,
				  const Environment& environment);

  double log_density;
  double log_density_rate;
  double density;
  double seeds_survival_weighted;
  double seeds_survival_weighted_rate;

  double pr_patch_survival_at_birth;

  static const int ode_dimension = 4;
};

}

RCPP_EXPOSED_CLASS_NODECL(model::CohortTop)

#endif
