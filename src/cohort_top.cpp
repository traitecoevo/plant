#include "cohort_top.h"
#include "gradient.h"

namespace model {

CohortTop::CohortTop(Strategy s) :
  Plant(s),
  density(0),
  density_rate(0),
  seeds_survival_weighted(0),
  seeds_survival_weighted_rate(0) {
}

CohortTop::CohortTop(Strategy *s) :
  Plant(s),
  density(0),
  density_rate(0),
  seeds_survival_weighted(0),
  seeds_survival_weighted_rate(0) {
}

// TODO: See design.md (search: compute_vars_phys_surv) for the
// issue around the name here.
void CohortTop::compute_vars_phys_surv(spline::Spline *env,
				       double survival_patch) {
  compute_vars_phys(env);

  // Defined on p. 7 at the moment.
  density_rate = growth_rate_gradient(env) + mortality_rate();

  // Defined on p 7, too.
  seeds_survival_weighted_rate =
    fecundity_rate() * survival_probability() * survival_patch;
}

// NOTE: germination_probability() will cause all physiological
// variables to be updated, so mass_leaf_rate() becomes valid so long
// as it is used afterwards.  This is something that can be improved,
// as that feels a bit fragile.
//
// NOTE: There will be a discussion of why the mortality rate initial
// condition is -log(germination_probability) in the documentation
// that Daniel is working out.
//
// NOTE: The initial condition for density is also a bit tricky, and
// defined on p 7 at the moment.
void CohortTop::compute_initial_conditions(spline::Spline *env,
					   double seed_input) {
  set_mortality(-log(germination_probability(env)));
  const double g = mass_leaf_rate();
  density = g > 0 ? (seed_input / g) : 0.0;
}

// * ODE interface
size_t CohortTop::ode_size() const {
  return ode_dimension;
}

ode::iter_const CohortTop::ode_values_set(ode::iter_const it) {
  set_mass_leaf(*it++);
  set_mortality(*it++);
  seeds_survival_weighted = *it++;
  density = *it++;
  return it;
}

ode::iter CohortTop::ode_values(ode::iter it) const {
  *it++ = mass_leaf();
  *it++ = mortality();
  *it++ = seeds_survival_weighted;
  *it++ = density;
  return it;
}

ode::iter CohortTop::ode_rates(ode::iter it) const {
  *it++ = mass_leaf_rate();
  *it++ = mortality_rate();
  *it++ = seeds_survival_weighted_rate;
  *it++ = density_rate;
  return it;
}

// * R interface
void CohortTop::r_compute_vars_phys_surv(spline::Spline env,
					 double survival_patch) {
  compute_vars_phys_surv(&env, survival_patch);
}
void CohortTop::r_compute_initial_conditions(spline::Spline env,
					     double seed_input) {
  compute_initial_conditions(&env, seed_input);
}
double CohortTop::r_growth_rate_gradient(spline::Spline env) const {
  return growth_rate_gradient(&env);
}
double CohortTop::r_growth_rate_given_mass(double mass_leaf,
					   spline::Spline env) {
  return growth_rate_given_mass(mass_leaf, &env);
}

// This is the gradient with respect to mass_leaf.  It is needed for
// computing the derivative (wrt time) of the density of individuals.
double CohortTop::growth_rate_gradient(spline::Spline *env) const {
  CohortTop tmp = *this;
  FunctorBind2<CohortTop, spline::Spline*,
	       &CohortTop::growth_rate_given_mass> fun(&tmp, env);

  const double d_mass_leaf = 1e-6; // TODO: Constant/untunable for now

  // const int r = 4;
  // return util::gradient_richardson(&fun, mass_leaf(), d_mass_leaf, r);

  return util::gradient_fd_forward(&fun, mass_leaf(),
				   d_mass_leaf, mass_leaf_rate());
}

// This exists only because it is needed by growth_rate_gradient.
double CohortTop::growth_rate_given_mass(double mass_leaf,
					 spline::Spline *env) {
  set_mass_leaf(mass_leaf);
  compute_vars_phys(env);
  return mass_leaf_rate();
}

}
