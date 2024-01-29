// -*-c++-*-
#ifndef PLANT_PLANT_K93_ENVIRONMENT_H_
#define PLANT_PLANT_K93_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>

using namespace Rcpp;

namespace plant {

class K93_Environment : public Environment {
public:
  K93_Environment() {
    time = 0.0;
    shading = Resource_spline();
  };

  // Light interface
  Resource_spline shading;

  // todo: Should this be here or in shading?
  void set_fixed_environment(double value, double height_max) {
    shading.set_fixed_values(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return shading.get_value_at_height(height);
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    shading.r_init_interpolators(state);
  }

  // Core functions
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of environment
    auto f_canopy_openness = [&](double height) -> double
    { return exp(-f_compute_competition(height)); };

    // Calculates the shading environment, fitting a spline to the the function
    // `f_compute_competition` as a function of height
    shading.compute_environment(f_canopy_openness, height_max, rescale);
  }

  virtual void clear_environment() {
    shading.clear();
  }

};


inline Rcpp::List get_state(const K93_Environment environment, double time) {
  auto ret = get_state(environment.extrinsic_drivers, time);
  ret["shading"] = get_state(environment.shading); // does a full copy of ret, not efficient
  return ret;
}

}

#endif
