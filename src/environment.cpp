#include "environment.h"

namespace model {

Environment::Environment(Parameters p)
  : disturbance_regime(p.mean_disturbance_interval),
    seed_rain(p.seed_rain),
    seed_rain_index(0),
    control(p.control),
    time(0.0),
    light_environment_generator(control.environment_light_tol,
				control.environment_light_tol,
				control.environment_light_nbase,
				control.environment_light_max_depth,
				control.environment_light_akima) {
}

double Environment::canopy_openness(double height) const {
  const bool within_canopy = height <= light_environment.max();
  return within_canopy ? light_environment.eval(height) : 1.0;
}

void Environment::compute_light_environment(util::DFunctor *f_canopy_openness,
					    double height_max) {
  light_environment =
    light_environment_generator.construct(f_canopy_openness,
					  0, height_max);
}

void Environment::rescale_light_environment(util::DFunctor *f_canopy_openness,
					    double height_max) {
  std::vector<double> h = light_environment.get_x();
  const double min = light_environment.min(),
    max_old = light_environment.max();
  light_environment.clear();

  util::rescale(h.begin(), h.end(), min, max_old, min, height_max);
  h[h.size() - 1] = height_max; // avoid round-off issues

  for (std::vector<double>::const_iterator hi = h.begin();
       hi != h.end(); ++hi) {
    light_environment.add_point(*hi, (*f_canopy_openness)(*hi));
  }

  light_environment.initialise();
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
  light_environment.clear();
}

double Environment::seed_rain_rate() const {
  if (seed_rain.empty())
    Rcpp::stop("Cannot get seed rain for empty environment");
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

const Disturbance& Environment::get_disturbance_regime() const {
  return disturbance_regime;
}

// * R interface
interpolator::Interpolator Environment::r_get_light_environment() const {
  return light_environment;
}
void Environment::r_set_light_environment(const interpolator::Interpolator env) {
  light_environment = env;
}

std::vector<double> Environment::r_get_seed_rain() const {
  return seed_rain;
}

void Environment::r_set_seed_rain_index(size_t x) {
  set_seed_rain_index(util::check_bounds_r(x, seed_rain.size()));
}

// The only *variable* in Enironment is the time.
//
// - seed_rain and disturbance are set by Parameters
// - light_environment follows from the patch
//
// However, because of control.rescale_light_environment, we might not
// always get the same light environment.
Rcpp::List Environment::r_get_state() const {
  return Rcpp::List::create(Rcpp::_["time"] = time);
}

void Environment::r_set_state(Rcpp::List x) {
  time = Rcpp::as<double>(x["time"]);
}


}
