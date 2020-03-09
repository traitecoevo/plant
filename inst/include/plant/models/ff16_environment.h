// -*-c++-*-
#ifndef PLANT_PLANT_LIGHT_ENVIRONMENT_H_
#define PLANT_PLANT_LIGHT_ENVIRONMENT_H_

#include <plant/control.h>
#include <plant/disturbance.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/environment.h>
#include <plant/util.h>

using namespace Rcpp;

namespace plant {

class FF16_Environment : public Environment {
public:

  FF16_Environment() {
    time = 0.0;
    disturbance_regime = 30;
    seed_rain = { 1.0, 1.0, 1.0 };
    seed_rain_index = 3;
    k_I = 0.5;
    environment_generator = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
    // Define an anonymous function to pass got the environment generator
    environment_interpolator = environment_generator.construct(
      [&](double height) {
          return get_environment_at_height(height);
      }, 0, 1); // these are update with init(x, y) whne patch is created
  };

  FF16_Environment(double disturbance_mean_interval,
                   std::vector<double> seed_rain_,
                   double k_I_,
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
  };

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
};

}

#endif
