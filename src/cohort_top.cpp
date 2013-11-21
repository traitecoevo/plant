#include "cohort_top.h"
#include "gradient.h"
#include "functor.h"

namespace model {

CohortTop::CohortTop(Strategy s) :
  Plant(s),
  log_density(R_NegInf),
  log_density_rate(0),
  density(0),
  seeds_survival_weighted(0),
  seeds_survival_weighted_rate(0),
  pr_patch_survival_at_birth(1) {
}

CohortTop::CohortTop(Strategy *s) :
  Plant(s),
  log_density(R_NegInf),
  log_density_rate(0),
  density(0),
  seeds_survival_weighted(0),
  seeds_survival_weighted_rate(0),
  pr_patch_survival_at_birth(1) {
}

void CohortTop::compute_vars_phys(const Environment& environment) {
  Plant::compute_vars_phys(environment);

  // [eqn 22] Density per unit area of individuals with size 'h'.
  // EBT.md{eq:boundN}, see Numerical technique.
  // details.md, for details on translation from mass_leaf to height.
  log_density_rate =
    - growth_rate_gradient(environment, get_last_integration_intervals())
    - mortality_rate();

  // EBT.md{eq:boundSurv}, see Numerical technique
  const double survival_patch = environment.patch_survival();
  seeds_survival_weighted_rate =
    fecundity_rate() * survival_probability() *
    survival_patch / pr_patch_survival_at_birth;
}

double CohortTop::leaf_area_above(double z) const {
  // EBT.md{eq:boundN}, and following section.
  return density * Plant::leaf_area_above(z);
}

double CohortTop::leaf_area() const {
  return density * Plant::leaf_area();
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
  pr_patch_survival_at_birth = environment.patch_survival();
  const double pr_germ = germination_probability(environment);
  // EBT.md{eq:boundSurv}
  set_mortality(-log(pr_germ));
  // EBT.md{eq:boundN}
  const double g = height_rate();
  const double seed_rain = environment.seed_rain_rate();
  log_density = log(g > 0 ? seed_rain * pr_germ / g : 0.0);
  density     = exp(log_density);
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
  density = exp(log_density);
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

double CohortTop::r_growth_rate_gradient(const Environment& environment) {
  Plant::compute_vars_phys(environment);
  return growth_rate_gradient(environment, get_last_integration_intervals());
}
double
CohortTop::r_growth_rate_given_height(double height_,
				      const Environment& environment) {
  return growth_rate_given_height(height_, environment);
}

size_t CohortTop::state_size() const {
  return ode_size() + 1;
}

CohortTop::state::iterator
CohortTop::get_state(CohortTop::state::iterator it) const {
  it = ode_values(it);
  *it++ = pr_patch_survival_at_birth;
  return it;
}

CohortTop::state::const_iterator
CohortTop::set_state(CohortTop::state::const_iterator it) {
  it = set_ode_values(0 /* unused - time */, it);
  pr_patch_survival_at_birth = *it++;
  return it;
}

// This is the gradient of height_rate with respect to height.  It is
// needed for computing the derivative (wrt time) of the log_density of
// individuals.
double CohortTop::growth_rate_gradient(const Environment& environment,
				       integration::intervals_type
				       intervals) const {
  CohortTop tmp = *this;
  tmp.set_integration_intervals(intervals);
  util::FunctorBind2<CohortTop, const Environment&,
		     &CohortTop::growth_rate_given_height>
    fun(&tmp, environment);

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
