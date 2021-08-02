// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/canopy.h>

using namespace Rcpp;

namespace plant {

class FF16_Environment : public Environment {
public:

  FF16_Environment();

  // Light interface
  double canopy_light_tol;
  size_t canopy_light_nbase;
  size_t canopy_light_max_depth;
  bool canopy_rescale_usually;

  Canopy canopy;

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

  // Should this be here or in canopy.h?
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


  // Soil interface
  double soil_infiltration_rate;
  size_t soil_number_of_depths;
  double soil_initial_state;

  Internals vars;

  // This is temporary
  virtual void compute_rates() {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.set_rate(i, soil_infiltration_rate / (i+1));
    }
  }

  std::vector<double> get_soil_water_state() const {
    return vars.states;
  }

  // R interface
  // This is temporary, eventually want to set vectors
  void set_soil_water_state(double v) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.set_state(i, v / (i+1));
    }
  }

  Internals r_internals() const { return vars; }
};

inline Rcpp::NumericMatrix get_state(const FF16_Environment environment) {
  return get_state(environment.canopy);
}


}

#endif
