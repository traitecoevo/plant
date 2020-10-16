// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_

#include <plant/control.h>
#include <plant/disturbance.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/ode_interface.h>
#include <plant/util.h>
#include <Rcpp.h>

using namespace Rcpp;

namespace plant {

class Environment {
public:

  template <typename Function>
  void compute_environment(Function f, double height_max);
  template <typename Function>
  void rescale_environment(Function f, double height_max);

  void set_fixed_environment(double value, double height_max);
  void set_fixed_environment(double value);
  void init_interpolators(const std::vector<double>& state);

  void compute_rates();

  // Dummy iterators: do nothing if the environment has no state. 
  ode::iterator ode_state(ode::iterator it) const { return it; }
  ode::const_iterator set_ode_state(ode::const_iterator it) { return it; }
  ode::iterator ode_rates(ode::iterator it) const { return it; }
  int ode_size() const { return 0; }

  // Computes the probability of survival from 0 to time.
  double patch_survival() const {
    return disturbance_regime.pr_survival(time);
  }

  // Computes the probability of survival from time_at_birth to time, by
  // conditioning survival over [0,time] on survival over
  // [0,time_at_birth].
  double patch_survival_conditional(double time_at_birth) const {
    return disturbance_regime.pr_survival_conditional(time, time_at_birth);
  }

  // Reset the environment.
  void clear() {
    time = 0.0;
    clear_environment();
  }

  void clear_environment() {}

  double seed_rain_dt() const {
    if (seed_rain.empty()) {
      Rcpp::stop("Cannot get seed rain for empty environment");
    }
    return seed_rain[seed_rain_index];
  }

  void set_seed_rain_index(size_t x) {
    seed_rain_index = x;
  }

  // * R interface
  void r_set_seed_rain_index(util::index x) {
    set_seed_rain_index(x.check_bounds(seed_rain.size()));
  }

  void r_init_interpolators(const std::vector<double>& state) {}

  double get_environment_at_height(double height) { return 0.0; };

  double time;
  Disturbance disturbance_regime;
  std::vector<double> seed_rain;
  size_t seed_rain_index;
};

}

#endif
