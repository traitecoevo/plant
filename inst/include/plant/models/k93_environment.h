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
    canopy = Canopy();
  };

  // Light interface
  bool shading_spline_rescale_usually;

  Canopy canopy;

  // Should this be here or in canopy?
  void set_fixed_environment(double value, double height_max) {
    canopy.set_fixed_canopy(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return canopy.get_canopy_at_height(height);
  }

  double canopy_openness(double height) const {
    return canopy.canopy_openness(height);
  }

  void r_init_interpolators(const std::vector<double>& state) {
    canopy.r_init_interpolators(state);
  }

  // Core functions
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max) {
    canopy.compute_canopy(f_compute_competition, height_max);
  }

  template <typename Function>
  void rescale_environment(Function f_compute_competition, double height_max) {
    canopy.rescale_canopy(f_compute_competition, height_max);
  }

  void clear_environment() {
    canopy.clear();
  }

};

//inline Rcpp::NumericMatrix get_state(const K93_Environment environment) {
//  return get_state(environment.canopy);
//}
inline Rcpp::List get_state(const K93_Environment environment, double time) {
  auto ret = get_state(environment.extrinsic_drivers, time);
  ret["canopy"] = get_state(environment.canopy); // does a full copy of ret, not efficient
  return ret;
}

}

#endif
