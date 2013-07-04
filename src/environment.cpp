#include "environment.h"

namespace model {

Environment::Environment(Parameters p)
  : disturbance_regime(p.disturbance_regime),
    seed_rain(p.size(), 0.0),
    seed_rain_index(0),
    control(p.control),
    time(0.0),
    light_environment_generator(control.environment_light_tol,
				control.environment_light_tol,
				control.environment_light_nbase,
				control.environment_light_max_depth) {
}

double Environment::canopy_openness(double height) const {
  const bool within_canopy = height <= light_environment.max();
  return within_canopy ? light_environment.eval(height) : 1.0;
}

void Environment::compute_light_environment(util::DFunctor *f_canopy_openness,
					    double height_max) {
  light_environment =
    light_environment_generator.construct_spline(f_canopy_openness,
						 0, height_max);
}

void Environment::scale_light_environment(util::DFunctor *f_canopy_openness,
					  double height_max) {
  std::vector<double> x = light_environment.get_x();
  std::vector<double> y;
  y.reserve(x.size());
  const double r = height_max / x.back();
  for (std::vector<double>::iterator xi = x.begin();
       xi != x.end(); ++xi) {
    *xi *= r;
    y.push_back((*f_canopy_openness)(*xi));
  }
  if (!util::identical(x.back(), height_max)) {
    ::Rf_warning("Can't rely on this it seems");
    x.back() = height_max;
    y.back() = (*f_canopy_openness)(height_max);
  }
  light_environment.init(x, y);
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
    ::Rf_error("Cannot get seed rain for empty environment");
  return seed_rain[seed_rain_index];
}

double Environment::get_time() const {
  return time;
}
void Environment::set_time(double x) {
  time = x;
}

void Environment::set_seed_rain_index(size_t x) {
  seed_rain_index = x;
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

void Environment::r_set_seed_rain_index(size_t x) {
  set_seed_rain_index(util::check_bounds_r(x, seed_rain.size()));
}


}
