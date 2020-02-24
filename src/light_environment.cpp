#include <plant/light_environment.h>

namespace plant {

LightEnvironment::LightEnvironment() {
  time = 0.0;
  disturbance_regime = 30;
  seed_rain = { 1.0, 1.0, 1.0 };
  seed_rain_index = 3;
  environment_generator = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
  // Define an anonymous function to pass got the environment generator
  environment_interpolator = environment_generator.construct(
    [&](double height) {
        return LightEnvironment::canopy_openness(height);
    }, 30, 150.0);
};


double LightEnvironment::canopy_openness(double height) const {
  const bool within_canopy = height <= environment_interpolator.max();
  return within_canopy ? environment_interpolator.eval(height) : 1.0;
}

// Computes the probability of survival from 0 to time.
double LightEnvironment::patch_survival() const {
  return disturbance_regime.pr_survival(time);
}

// Computes the probability of survival from time_at_birth to time, by
// conditioning survival over [0,time] on survival over
// [0,time_at_birth].
double LightEnvironment::patch_survival_conditional(double time_at_birth) const {
  return disturbance_regime.pr_survival_conditional(time, time_at_birth);
}

// Reset the environment.
void LightEnvironment::clear() {
  time = 0.0;
  clear_environment();
}

void LightEnvironment::clear_environment() {
  environment_interpolator.clear();
}

double LightEnvironment::seed_rain_dt() const {
  if (seed_rain.empty()) {
    Rcpp::stop("Cannot get seed rain for empty environment");
  }
  return seed_rain[seed_rain_index];
}

void LightEnvironment::set_seed_rain_index(size_t x) {
  seed_rain_index = x;
}

// * R interface
void LightEnvironment::r_set_seed_rain_index(util::index x) {
  set_seed_rain_index(x.check_bounds(seed_rain.size()));
}

}
