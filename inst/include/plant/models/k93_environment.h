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
  };

  // Light interface
  ResourceSpline<double> light_availability;

  void set_fixed_environment(double value, double height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return light_availability.get_value_at_height(height);
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }
  virtual std::vector<double> get_interpolators_state() const {
    return light_availability.interpolators_state();
  }

  virtual Rcpp::List r_get_state() const
  {
    return Rcpp::List::create(_["light_availability"] = time); //      light_availability);
  }

  // Core functions
  // The light a height is left with, from the competition profile above it.
  template <typename Function>
  void compute_environment(Function f_compute_competition_and_slope,
                           double height_max) {
    build_extinction_field(light_availability, f_compute_competition_and_slope,
                           height_max);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }

};

}

#endif
