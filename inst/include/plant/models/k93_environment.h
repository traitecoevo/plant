// -*-c++-*-
#ifndef PLANT_PLANT_K93_ENVIRONMENT_H_
#define PLANT_PLANT_K93_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/shading.h>

using namespace Rcpp;

namespace plant {

class K93_Environment : public Environment {
public:
  K93_Environment() {
    shading_spline_rescale_usually = false;
    time = 0.0;
    shading = Shading();
  };

  // Light interface
  bool shading_spline_rescale_usually;

  Shading shading;

  // todo: Should this be here or in shading?
  void set_fixed_environment(double value, double height_max) {
    shading.set_fixed_values(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return shading.get_canopy_at_height(height);
  }

  double canopy_openness(double height) const {
    return shading.canopy_openness(height);
  }

  void r_init_interpolators(const std::vector<double>& state) {
    shading.r_init_interpolators(state);
  }

  // Core functions
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max) {
    shading.compute_shading(f_compute_competition, height_max);
  }

  template <typename Function>
  void rescale_environment(Function f_compute_competition, double height_max) {
    shading.rescale_points(f_compute_competition, height_max);
  }

  void clear_environment() {
    shading.clear();
  }

};

//inline Rcpp::NumericMatrix get_state(const K93_Environment environment) {
//  return get_state(environment.shading);
//}
inline Rcpp::List get_state(const K93_Environment environment, double time) {
  auto ret = get_state(environment.extrinsic_drivers, time);
  ret["shading"] = get_state(environment.shading); // does a full copy of ret, not efficient
  return ret;
}

}

#endif
