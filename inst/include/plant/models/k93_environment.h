// Built from  inst/include/plant/models/ff16r_environment.h on Wed Aug 12 15:46:38 2020 using the scaffolder, from the strategy:  FF16r
// Built from  inst/include/plant/models/ff16_environment.h on Wed Aug 12 11:12:34 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_K93_ENVIRONMENT_H_
#define PLANT_PLANT_K93_ENVIRONMENT_H_

#include <plant/control.h>
#include <plant/disturbance.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/environment.h>
#include <plant/util.h>
#include <Rcpp.h>

using namespace Rcpp;

namespace plant {

class K93_Environment : public Environment {
public:

  K93_Environment() {
    // Define an anonymous function to pass got the environment generator
    time = NA_REAL;
    disturbance_regime = 0;
    seed_rain = { 1.0, 1.0, 1.0 };
    seed_rain_index = 0;
    k_I = NA_REAL;
    environment_generator = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
    environment_interpolator = environment_generator.construct(
      [&](double height) {
          return get_environment_at_height(height);
      }, 0, 1); // these are update with init(x, y) when patch is created
  };

  K93_Environment(double disturbance_mean_interval,
                   std::vector<double> seed_rain_,
                   double k_I_,
                   Control control) {
    k_I = k_I_;
    time = 0.0;
    disturbance_regime = disturbance_mean_interval;
    seed_rain = seed_rain_;
    seed_rain_index = 0;
    environment_generator = interpolator::AdaptiveInterpolator(
      control.environment_light_tol,
      control.environment_light_tol,
      control.environment_light_nbase,
      control.environment_light_max_depth
    );
    environment_interpolator = environment_generator.construct(
      [&](double height) {
          return get_environment_at_height(height);
      }, 0, 1); // these are update with init(x, y) when patch is created
  };

  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max) {
    const double lower_bound = 0.0;
    double upper_bound = height_max;

    auto f_canopy_openness = [&] (double height) -> double {return exp(-k_I * f_compute_competition(height));};
    environment_interpolator =
      environment_generator.construct(f_canopy_openness, lower_bound, upper_bound);
  }

  template <typename Function>
  void rescale_environment(Function f_compute_competition, double height_max) {
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

  void set_fixed_environment(double value, double height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<double> y = {value, value, value};
    clear_environment();
    environment_interpolator.init(x, y);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double canopy_openness(double height) const {
    return get_environment_at_height(height);
  }
};

inline Rcpp::NumericMatrix get_state(const K93_Environment environment) {
  using namespace Rcpp;
  NumericMatrix xy = environment.environment_interpolator.r_get_xy();
  Rcpp::CharacterVector colnames =
    Rcpp::CharacterVector::create("height", "canopy_openness");
  xy.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);
  return xy;
}


}

#endif
