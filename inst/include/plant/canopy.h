// -*-c++-*-
#ifndef PLANT_PLANT_CANOPY_H_
#define PLANT_PLANT_CANOPY_H_

#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/ode_interface.h>
#include <plant/util.h>

using namespace Rcpp;

namespace plant {

class Canopy {
public:

  Canopy() {
    // Define an anonymous function to pass got the environment generator
    k_I = NA_REAL;
    canopy_generator = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
    canopy_interpolator = canopy_generator.construct(
      [&](double height) {
          return get_canopy_at_height(height);
      }, 0, 1); // these are update with init(x, y) when patch is created
  };

  Canopy(double k_I_, Control control) {
    k_I = k_I_;
    canopy_generator = interpolator::AdaptiveInterpolator(
      control.environment_light_tol,
      control.environment_light_tol,
      control.environment_light_nbase,
      control.environment_light_max_depth
    );
    canopy_interpolator = canopy_generator.construct(
      [&](double height) {
          return get_canopy_at_height(height);
      }, 0, 1); // these are update with init(x, y) when patch is created
  };

  template <typename Function>
  void compute_canopy(Function f_compute_competition, double height_max) {
    const double lower_bound = 0.0;
    double upper_bound = height_max;

    auto f_canopy_openness = [&] (double height) -> double {return exp(-k_I * f_compute_competition(height));};
    canopy_interpolator =
      canopy_generator.construct(f_canopy_openness, lower_bound, upper_bound);
  }

  template <typename Function>
  void rescale_canopy(Function f_compute_competition, double height_max) {
    std::vector<double> h = canopy_interpolator.get_x();
    const double min = canopy_interpolator.min(), // 0.0?
      height_max_old = canopy_interpolator.max();

    auto f_canopy_openness = [&] (double height) -> double {return exp(-k_I * f_compute_competition(height));};
    util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
    h.back() = height_max; // Avoid round-off error.

    canopy_interpolator.clear();
    for (auto hi : h) {
      canopy_interpolator.add_point(hi, f_canopy_openness(hi));
    }
    canopy_interpolator.initialise();
  }

  void set_fixed_canopy(double value, double height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<double> y = {value, value, value};
    clear();
    canopy_interpolator.init(x, y);
  }

  void clear() {
    canopy_interpolator.clear();
  }

  double get_canopy_at_height(double height) const {
    const bool within = height <= canopy_interpolator.max();
    // TODO: change maximum - here hard-coded to 1.0
    return within ? canopy_interpolator.eval(height) : 1.0;
  }

  double canopy_openness(double height) const {
    return get_canopy_at_height(height);
  }


  void r_init_interpolators(const std::vector<double>& state) {
    // See issue #144; this is important as we have to at least refine
    // the light environment, but doing this is better because it means
    // that if rescale_usually is on we do get the same light
    // environment as before.
    if (state.size() % 2 != 0) {
      util::stop("Expected even number of elements in light environment");
    }
    const size_t state_n = state.size() / 2;
    auto it = state.begin();
    std::vector<double> state_x, state_y;
    std::copy_n(it,         state_n, std::back_inserter(state_x));
    std::copy_n(it + state_n, state_n, std::back_inserter(state_y));
    canopy_interpolator.init(state_x, state_y);
  }

  double k_I;
  interpolator::Interpolator canopy_interpolator;
  interpolator::AdaptiveInterpolator canopy_generator;
};

inline Rcpp::NumericMatrix get_state(const Canopy canopy) {
  using namespace Rcpp;
  NumericMatrix xy = canopy.canopy_interpolator.r_get_xy();
  Rcpp::CharacterVector colnames =
    Rcpp::CharacterVector::create("height", "canopy_openness");
  xy.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);
  return xy;
}


}

#endif
