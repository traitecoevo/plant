#include <plant/light_environment.h>

namespace plant {

FF16_Environment::FF16_Environment() {
  time = 0.0;
  disturbance_regime = 30;
  seed_rain = { 1.0, 1.0, 1.0 };
  seed_rain_index = 3;
  environment_generator = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
  // Define an anonymous function to pass got the environment generator
  environment_interpolator = environment_generator.construct(
    [&](double height) {
        return FF16_Environment::canopy_openness(height);
    }, 30, 150.0);
  k_I = 0.5;
};

void FF16_Environment::set_fixed_environment(double canopy_openness, double height_max) {
  std::vector<double> x = {0, height_max/2.0, height_max};
  std::vector<double> y = {canopy_openness, canopy_openness, canopy_openness};
  clear_environment();
  environment_interpolator.init(x, y);
}

void FF16_Environment::set_fixed_environment(double canopy_openness) {
  double height_max = 150.0;
  set_fixed_environment(canopy_openness, height_max);
}


double FF16_Environment::canopy_openness(double height) const {
  const bool within_canopy = height <= environment_interpolator.max();
  return within_canopy ? environment_interpolator.eval(height) : 1.0;
}

// Computes the probability of survival from 0 to time.
double FF16_Environment::patch_survival() const {
  return disturbance_regime.pr_survival(time);
}

// Computes the probability of survival from time_at_birth to time, by
// conditioning survival over [0,time] on survival over
// [0,time_at_birth].
double FF16_Environment::patch_survival_conditional(double time_at_birth) const {
  return disturbance_regime.pr_survival_conditional(time, time_at_birth);
}

// Reset the environment.
void FF16_Environment::clear() {
  time = 0.0;
  clear_environment();
}

void FF16_Environment::clear_environment() {
  environment_interpolator.clear();
}

double FF16_Environment::seed_rain_dt() const {
  if (seed_rain.empty()) {
    Rcpp::stop("Cannot get seed rain for empty environment");
  }
  return seed_rain[seed_rain_index];
}

void FF16_Environment::set_seed_rain_index(size_t x) {
  seed_rain_index = x;
}

// * R interface
void FF16_Environment::r_set_seed_rain_index(util::index x) {
  set_seed_rain_index(x.check_bounds(seed_rain.size()));
}

}
