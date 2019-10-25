#include <plant/environment.h>
#include <plant/parameters.h>

namespace plant {

Environment::Environment(double disturbance_mean_interval,
                         std::vector<double> seed_rain_,
                         Control control)
  : time(0.0),
    disturbance_regime(disturbance_mean_interval),
    seed_rain(seed_rain_),
    seed_rain_index(0),
    light_environment_generator(make_interpolator(control)) {
}

double Environment::canopy_openness(double size) const {
  const bool within_canopy = size <= light_environment.max();
  return within_canopy ? light_environment.eval(size) : 1.0;
}

// Computes the probability of survival from 0 to time.
double Environment::patch_survival() const {
  return disturbance_regime.pr_survival(time);
}

// Computes the probability of survival from time_at_birth to time, by
// conditioning survival over [0,time] on survival over
// [0,time_at_birth].
double Environment::patch_survival_conditional(double time_at_birth) const {
  return disturbance_regime.pr_survival_conditional(time, time_at_birth);
}

// Reset the environment.
void Environment::clear() {
  time = 0.0;
  clear_light_environment();
}

void Environment::clear_light_environment() {
  light_environment.clear();
}

double Environment::seed_rain_dt() const {
  if (seed_rain.empty()) {
    Rcpp::stop("Cannot get seed rain for empty environment");
  }
  return seed_rain[seed_rain_index];
}

void Environment::set_seed_rain_index(size_t x) {
  seed_rain_index = x;
}

// * R interface
void Environment::r_set_seed_rain_index(util::index x) {
  set_seed_rain_index(x.check_bounds(seed_rain.size()));
}

}
