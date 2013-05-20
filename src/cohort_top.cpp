#include "cohort_top.h"
#include "gradient.h"

namespace model {

CohortTop::CohortTop(Strategy s) :
  Plant(s),
  density(0),
  density_rate(0),
  seeds_survival_weighted(0),
  seeds_survival_weighted_rate(0),
  time_of_birth(0) {
}

CohortTop::CohortTop(Strategy *s) :
  Plant(s),
  density(0),
  density_rate(0),
  seeds_survival_weighted(0),
  seeds_survival_weighted_rate(0),
  time_of_birth(0) {
}

// TODO: See design.md (search: compute_vars_phys_surv) for the
// issue around the name here.
void CohortTop::compute_vars_phys(const Environment& environment) {
  Plant::compute_vars_phys(environment);

  // Defined on p. 7 at the moment.
  density_rate = growth_rate_gradient(environment) + mortality_rate();

  // Defined on p 7, too.
  const double survival_patch = environment.patch_survival(time_of_birth);
  seeds_survival_weighted_rate =
    fecundity_rate() * survival_probability() * survival_patch;
}

// TODO: Not sure if this should happen here, of we should treat this
// as the *individual* leaf area above?
double CohortTop::leaf_area_above(double z) const {
  return density * Plant::leaf_area_above(z);
}

int CohortTop::offspring() {
  ::Rf_error("Cannot use offspring() in CohortTop");
  return 0;
}

bool CohortTop::died() {
  ::Rf_error("Cannot use died() in CohortTop");
  return false;
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
void CohortTop::compute_initial_conditions(const Environment& environment) {
  time_of_birth = environment.get_age();
  set_mortality(-log(germination_probability(environment)));
  const double g = mass_leaf_rate();
  const double seed_rain = environment.seed_rain_rate();
  density = g > 0 ? (seed_rain / g) : 0.0;
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

double CohortTop::r_growth_rate_gradient(const Environment& environment) 
  const {
  return growth_rate_gradient(environment);
}
double CohortTop::r_growth_rate_given_mass(double mass_leaf, 
					   const Environment& environment) {
  return growth_rate_given_mass(mass_leaf, environment);
}

// This is the gradient of mass_leaf_rate with respect to mass_leaf.
// It is needed for computing the derivative (wrt time) of the density
// of individuals.
double CohortTop::growth_rate_gradient(const Environment& environment) const {
  CohortTop tmp = *this;
  FunctorBind2<CohortTop, const Environment&,
	       &CohortTop::growth_rate_given_mass> fun(&tmp, environment);

  const double eps = control().cohort_gradient_eps;
  double grad;
  if (control().cohort_gradient_richardson) {
    const int r = control().cohort_gradient_richardson_depth;
    grad = util::gradient_richardson(&fun, mass_leaf(), eps, r);
  } else {
    grad = util::gradient_fd_forward(&fun, mass_leaf(),
				     eps, mass_leaf_rate());
  }
  return grad;
}

// This exists only because it is needed by growth_rate_gradient.
double CohortTop::growth_rate_given_mass(double mass_leaf,
					 const Environment& environment) {
  set_mass_leaf(mass_leaf);
  Plant::compute_vars_phys(environment);
  return mass_leaf_rate();
}

}
