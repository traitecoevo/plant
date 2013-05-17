#include "environment.h"

#include "adaptive_spline.h"

namespace model {

Environment::Environment(Parameters p)
  : disturbance_regime(p.disturbance_regime),
    control(p.control),
    age(0.0) {
}

double Environment::canopy_openness(double height) const {
  const bool within_canopy = height <= light_environment.max();
  return within_canopy ? light_environment.eval(height) : 1.0;
}

void Environment::compute_light_environment(util::DFunctor *canopy_openness,
					    double height_max) {
  spline::AdaptiveSpline generator(canopy_openness);
  generator.set_control(control.environment_light_tol,
			control.environment_light_tol,
			control.environment_light_nbase,
			control.environment_light_max_depth);
  set_light_environment(generator.construct_spline(0, height_max));
}

// Computes the probability of survival from age_at_birth to age,
// probably by conditioning survival over [0,age] on survival over
// [0,age_at_birth].
double Environment::patch_survival(double age_at_birth) const {
  return disturbance_regime.survival_probability(age_at_birth, age);
}

// TODO: If I'm just going to get and set, should these just be public
// data?  What is the advantage here?  Is anything actually being
// encapsulated?
spline::Spline Environment::get_light_environment() const {
  return light_environment;
}
void Environment::set_light_environment(const spline::Spline env) {
  light_environment = env;
}

double Environment::get_age() const {
  return age;
}
void Environment::set_age(double x) {
  age = x;
}

}
