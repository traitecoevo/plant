#include "cohort_top.h"
#include "gradient.h"

namespace model {

CohortTop::CohortTop(Strategy s) :
  Plant(s),
  log_density(R_NegInf),
  log_density_rate(0),
  seeds_survival_weighted(0),
  seeds_survival_weighted_rate(0),
  time_of_birth(0) {
}

CohortTop::CohortTop(Strategy *s) :
  Plant(s),
  log_density(R_NegInf),
  log_density_rate(0),
  seeds_survival_weighted(0),
  seeds_survival_weighted_rate(0),
  time_of_birth(0) {
}

void CohortTop::compute_vars_phys(const Environment& environment) {
  Plant::compute_vars_phys(environment);

  // [eqn 22] Density per unit area of individuals with size 'h'.
  // EBT.md{eq:boundN}, see Numerical technique.
  // details.md, for details on translation from mass_leaf to height.
  log_density_rate = -growth_rate_gradient(environment) - mortality_rate();

  // EBT.md{eq:boundSurv}, see Numerical technique
  const double survival_patch = environment.patch_survival(time_of_birth);
  seeds_survival_weighted_rate =
    fecundity_rate() * survival_probability() * survival_patch;
}

double CohortTop::leaf_area_above(double z) const {
  // EBT.md{eq:boundN}, and following section.
  return exp(log_density) * Plant::leaf_area_above(z);
}

int CohortTop::offspring() {
  ::Rf_error("Cannot use offspring() in CohortTop");
  // return 0; // unreachable because Rf_error is noreturn (3.0.1)
}

bool CohortTop::died() {
  ::Rf_error("Cannot use died() in CohortTop");
  // return false; // unreachable because Rf_error is noreturn (3.0.1)
}

// NOTE: germination_probability() will cause all physiological
// variables to be updated, so height_rate() becomes valid so long
// as it is used afterwards.  This is something that can be improved,
// as that feels a bit fragile.
//
// NOTE: There will be a discussion of why the mortality rate initial
// condition is -log(germination_probability) in the documentation
// that Daniel is working out.
//
// NOTE: The initial condition for log_density is also a bit tricky, and
// defined on p 7 at the moment.
void CohortTop::compute_initial_conditions(const Environment& environment) {
  time_of_birth = environment.get_time();
  const double pr_germ = germination_probability(environment);
  // EBT.md{eq:boundSurv}
  set_mortality(-log(pr_germ));
  // EBT.md{eq:boundN}
  const double g = height_rate();
  const double seed_rain = environment.seed_rain_rate();
  log_density = log(g > 0 ? seed_rain * pr_germ / g : 0.0);
}

// * ODE interface
size_t CohortTop::ode_size() const {
  return ode_dimension;
}

ode::iterator_const CohortTop::set_ode_values(double /* unused: time */,
					      ode::iterator_const it) {
  set_height(*it++);
  set_mortality(*it++);
  seeds_survival_weighted = *it++;
  log_density = *it++;
  return it;
}

ode::iterator CohortTop::ode_values(ode::iterator it) const {
  *it++ = height();
  *it++ = mortality();
  *it++ = seeds_survival_weighted;
  *it++ = log_density;
  return it;
}

ode::iterator CohortTop::ode_rates(ode::iterator it) const {
  *it++ = height_rate();
  *it++ = mortality_rate();
  *it++ = seeds_survival_weighted_rate;
  *it++ = log_density_rate;
  return it;
}

double CohortTop::r_growth_rate_gradient(const Environment& environment)
  const {
  return growth_rate_gradient(environment);
}

// This is the gradient of height_rate with respect to height.  It is
// needed for computing the derivative (wrt time) of the log_density of
// individuals.
double CohortTop::growth_rate_gradient(const Environment& environment) const {
  CohortTop tmp = *this;
  FunctorBind2<CohortTop, const Environment&,
	       &CohortTop::growth_rate_given_height> fun(&tmp, environment);

  const double eps = control().cohort_gradient_eps;
  double grad;
  if (control().cohort_gradient_richardson) {
    const size_t r = control().cohort_gradient_richardson_depth;
    grad = util::gradient_richardson(&fun, height(), eps, r);
  } else if (control().cohort_gradient_direction > 0) {
    grad = util::gradient_fd_forward(&fun, height(),
				     eps, height_rate());
  } else {
    grad = util::gradient_fd_backward(&fun, height(),
				      eps, height_rate());
  }
  return grad;
}

// This exists only because it is needed by growth_rate_gradient.
double CohortTop::growth_rate_given_height(double height_,
					   const Environment& environment) {
  set_height(height_);
  Plant::compute_vars_phys(environment);
  return height_rate();
}

}
