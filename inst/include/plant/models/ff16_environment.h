// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/canopy.h>

using namespace Rcpp;

namespace plant {

class FF16_Environment : public Environment {
public:

  FF16_Environment() {
    // Define an anonymous function to pass got the environment generator
    time = NA_REAL;
    birth_rate = { 1.0, 1.0, 1.0 };
    species_arriving_index = 0;
    canopy = Canopy();
  };

  FF16_Environment(double disturbance_mean_interval,
                   std::vector<double> birth_rate_,
                   Control control) {
    time = 0.0;
    birth_rate = birth_rate_;
    species_arriving_index = 0;
    canopy = Canopy(control);
  };

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

  void set_fixed_environment(double value, double height_max) {
    canopy.set_fixed_canopy(value, height_max);
  }

  // Should this be here or in canopy?
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

  Canopy canopy;
};

inline Rcpp::NumericMatrix get_state(const FF16_Environment environment) {
  return get_state(environment.canopy);
}


}

#endif
