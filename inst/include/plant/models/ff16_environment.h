// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/canopy.h>
#include <plant/interpolator.h>
#include <iostream>

using namespace Rcpp;

namespace plant {

class FF16_Environment : public Environment {
public:
  // constructor for R interface - default settings can be modified
  // except for soil_number_of_depths and canopy_rescale_usually
  // which are only updated on construction
  FF16_Environment(bool canopy_rescale_usually = false,
                   int soil_number_of_depths = 0)
      : canopy_rescale_usually(canopy_rescale_usually) {
    time = 0.0;
    canopy = Canopy();
    vars = Internals(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, 0.0));
    extrinsic_drivers["rainfall"] = interpolator::Interpolator();
  };


  // Light interface
  bool canopy_rescale_usually;

  // private?
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

  // currently only used by FF16w but present in both FF16 and K93
  // avoids uninitialised variables/splines because vars.state_size = 0
  // perhaps not completely user safe, but I think diffcult to break from R API
  virtual void compute_rates(std::vector<std::vector<double>> const& resource_depletion) {
    // treat each soil layer as a separate resource pool
    for (size_t i = 0; i < vars.state_size; i++) {
      double water_infiltration = extrinsic_drivers["rainfall"].eval(time) / (i+1);

      double water_exfiltration;

      // sum species to get the total depletion at each layer
      for(size_t j = 0; resource_depletion.size(); j++) {
        water_exfiltration += resource_depletion[j][i];
      }
      double water_balance = water_infiltration - water_exfiltration;

      vars.set_rate(i, water_balance);
    }
  }

  std::vector<double> get_soil_water_state() const {
    return vars.states;
  }

  // I wonder if this needs a better name? See also environment.h
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

inline Rcpp::NumericMatrix get_state(const FF16_Environment environment) {
  return get_state(environment.canopy);
}


}

#endif
