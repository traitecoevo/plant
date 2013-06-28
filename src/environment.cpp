#include "environment.h"

#include "adaptive_spline.h"

namespace model {

Environment::Environment(Parameters p)
  : disturbance_regime(p.disturbance_regime),
    seed_rain(p.size(), 0.0),
    control(p.control),
    time(0.0) {
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
  light_environment = generator.construct_spline(0, height_max);
}

// Computes the probability of survival from time_at_birth to time,
// probably by conditioning survival over [0,time] on survival over
// [0,time_at_birth].
double Environment::patch_survival(double time_at_birth) const {
  return disturbance_regime.survival_probability(time_at_birth, time);
}

// Reset the environment.
void Environment::clear() {
  time = 0.0;
  light_environment.clear();
}

double Environment::seed_rain_rate() const {
  if (seed_rain.size() == 0)
    ::Rf_error("Requested rain out of bounds");
  return seed_rain[0];
}

double Environment::get_time() const {
  return time;
}
void Environment::set_time(double x) {
  time = x;
}

// * R interface
spline::Spline Environment::r_get_light_environment() const {
  return light_environment;
}
void Environment::r_set_light_environment(const spline::Spline env) {
  light_environment = env;
}

std::vector<double> Environment::r_get_seed_rain() const {
  return seed_rain;
}

void Environment::r_set_seed_rain(std::vector<double> x) {
  util::check_length(x.size(), seed_rain.size());
  seed_rain = x;
}


}
