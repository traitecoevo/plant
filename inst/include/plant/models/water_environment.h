// -*-c++-*-
#ifndef PLANT_PLANT_WATER_ENVIRONMENT_H_
#define PLANT_PLANT_WATER_ENVIRONMENT_H_

#include <plant/environment.h>

using namespace Rcpp;

namespace plant {

class Water_Environment : public Environment {
public:

  const int SOIL_WATER_INDEX = 1;

  Water_Environment() {
    time = NA_REAL;
    disturbance_regime = 30;
    seed_rain = { 1.0, 1.0, 1.0 };
    seed_rain_index = 0;
    k_I = NA_REAL;
    environment_generator = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
    // Define an anonymous function to pass got the environment generator
    environment_interpolator = environment_generator.construct(
      [&](double height) {
          return get_environment_at_height(height);
      }, 0, 1); // these are update with init(x, y) whne patch is created
    inflow_rate = 1.0;
    vars = Internals(1);
  };

  Water_Environment(double disturbance_mean_interval,
                   std::vector<double> seed_rain_, double k_I_,
                   Control control) {
    k_I = k_I_;
    environment_generator = interpolator::AdaptiveInterpolator(control.environment_light_tol,
                              control.environment_light_tol,
                              control.environment_light_nbase,
                              control.environment_light_max_depth);
    time = 0.0;
    disturbance_regime = disturbance_mean_interval;
    seed_rain = seed_rain_;
    seed_rain_index = 0;
    inflow_rate = 1.0;
    vars = Internals(1);
  };

  static size_t state_size () { return 1; }
  static size_t ode_size() { return state_size(); }

  void compute_rates() {
    set_soil_water_state(inflow_rate);
  }

  void set_soil_water_state(double v) {
    vars.set_state(SOIL_WATER_INDEX, v);
  }

  double get_soil_water_state() {
    vars.state(SOIL_WATER_INDEX);
  }

  ode::const_iterator set_ode_state(ode::const_iterator it) {
    for (int i = 0; i < vars.state_size; i++) {
      vars.states[i] = *it++;
    }
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    for (int i = 0; i < vars.state_size; i++) {
      *it++ = vars.states[i];
    }
    return it;
  }
  ode::iterator ode_rates(ode::iterator it) const {
    for (int i = 0; i < vars.state_size; i++) {
      *it++ = vars.rates[i];
    }
    return it;
  }


  // TODO: move these to Environment 
  template <typename Function>
  void compute_environment(Function f_compute_competition,
                           double height_max) {
    const double lower_bound = 0.0;
    double upper_bound = height_max;

    auto f_canopy_openness = [&] (double height) -> double {return exp(-k_I * f_compute_competition(height));};
    environment_interpolator =
      environment_generator.construct(f_canopy_openness, lower_bound, upper_bound);
  }

  template <typename Function>
  void rescale_environment(Function f_compute_competition,
                                             double height_max) {
    std::vector<double> h = environment_interpolator.get_x();
    const double min = environment_interpolator.min(), // 0.0?
      height_max_old = environment_interpolator.max();

    auto f_canopy_openness = [&] (double height) -> double {return exp(-k_I * f_compute_competition(height));};
    util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
    h.back() = height_max; // Avoid round-off error.

    environment_interpolator.clear();
    for (auto hi : h) {
      environment_interpolator.add_point(hi, f_canopy_openness(hi));
    }
    environment_interpolator.initialise();
  }

  double canopy_openness(double height) const {
    return get_environment_at_height(height);
  }

  double k_I;
  double inflow_rate;

  Internals r_internals() const { return vars; }

  Internals vars;
};

}

#endif
