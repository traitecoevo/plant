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
    environment_generator = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
    // Define an anonymous function to pass got the environment generator
    environment_interpolator = environment_generator.construct(
      [&](double height) {
          return FF16_Environment::canopy_openness(height);
      }, 30, 150.0);
    k_I = 0.5;
  };

  FF16_Environment(double disturbance_mean_interval,
                   std::vector<double> seed_rain_,
                   double k_I_,
                   Control control)
    : time(0.0),
      disturbance_regime(disturbance_mean_interval),
      seed_rain(seed_rain_),
      seed_rain_index(0),
      k_I(k_I_),
      environment_generator(interpolator::AdaptiveInterpolator(control.environment_light_tol,
                                control.environment_light_tol,
                                control.environment_light_nbase,
                                control.environment_light_max_depth)) {
  };

  void set_fixed_environment(double competition_amount, double height_max);
  void set_fixed_environment(double competition_amount);
  double canopy_openness(double height) const;
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max);
  template <typename Function>
  void rescale_environment(Function f_compute_competition, double height_max);
  double patch_survival() const;
  double patch_survival_conditional(double time_at_birth) const;
  void clear();
  void clear_environment();

  // NOTE: Interface here will change
  double seed_rain_dt() const;
  void set_seed_rain_index(size_t x);
  // * R interface
  void r_set_seed_rain_index(util::index x);

  double time;
  Disturbance disturbance_regime;
  interpolator::Interpolator environment_interpolator;

private:
  std::vector<double> seed_rain;
  size_t seed_rain_index;
  interpolator::AdaptiveInterpolator environment_generator;
  double k_I;
};

template <typename Function>
void FF16_Environment::compute_environment(Function f_compute_competition,
                                            double height_max) {
  const double lower_bound = 0.0;
  double upper_bound = height_max;

  auto f_canopy_openness = [&] (double height) -> double {return exp(-k_I * f_compute_competition(height));};
  environment_interpolator =
    environment_generator.construct(f_canopy_openness, lower_bound, upper_bound);
}

template <typename Function>
void FF16_Environment::rescale_environment(Function f_compute_competition,
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

}

#endif
