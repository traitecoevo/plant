// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/shading.h>
#include <plant/interpolator.h>

using namespace Rcpp;

namespace plant {

class FF16_Environment : public Environment {
public:
  // constructor for R interface - default settings can be modified
  // except for soil_number_of_depths and shading_spline_rescale_usually
  // which are only updated on construction
  FF16_Environment(bool shading_spline_rescale_usually = false,
                   int soil_number_of_depths = 0)
      : shading_spline_rescale_usually(shading_spline_rescale_usually) {
    time = 0.0;
    shading = Shading();
    vars = Internals(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, 0.0));
  };


  // Light interface

  // Shading object is used for calculating shading
  Shading shading;
  bool shading_spline_rescale_usually;

  // Ability to prescribe a fixed value
  // TODO: add setting to set other variables like water
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

  // TODO: shading.XXX_openness = shading.get_canopy_at_height
  double XXX_openness(double height) const {
    return shading.XXX_openness(height);
  }

  void r_init_interpolators(const std::vector<double>& state) {
    shading.r_init_interpolators(state);
  }

  virtual void compute_rates(std::vector<double> const& resource_depletion) {
    double infiltration;
    double net_flux;

    double drainage_multiplier = 0.1; // experimental only;

    // treat each soil layer as a separate resource pool
    for (size_t i = 0; i < vars.state_size; i++) {

      // initial representation of drainage; to be improved
      if(i == 0) {
        infiltration = extrinsic_drivers.evaluate("rainfall", time);
      } else {
        infiltration = std::max(vars.state(i - 1), 0.0) * drainage_multiplier;
      }

      // ecologically, soil water shouldn't go below zero
      // truncating at zero until such a model is implemented
      double drainage_rate = std::max(vars.state(i), 0.0) * drainage_multiplier;

      net_flux = infiltration - resource_depletion[i] - drainage_rate;
      vars.set_rate(i, net_flux);
    }

  }

  std::vector<double> get_soil_water_state() const {
    return vars.states;
  }

  // TODO: I wonder if this needs a better name? See also environment.h
  Internals r_internals() const { return vars; }

  // R interface
  void set_soil_water_state(std::vector<double> state) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.set_state(i, state[i]);
    }
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

//inline Rcpp::NumericMatrix get_state(const FF16_Environment environment) {
//  return get_state(environment.shading);
//}
inline Rcpp::List get_state(const FF16_Environment environment, double time) {
  auto ret = get_state(environment.extrinsic_drivers, time);
  ret["shading"] = get_state(environment.shading); // does a full copy of ret, not efficient
  return ret;
}
}

#endif
