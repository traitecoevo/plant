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
    rainfall_spline = interpolator::Interpolator();
  };

  // first spline computation
  void rainfall_init(std::vector<double> const& x, std::vector<double> const& y) {
    rainfall_spline.init(x, y); // set x, y points and compute spline
  }

  // if any points are added to the spline, need to recompute spline
  void rainfall_recompute() {
    rainfall_spline.initialise();
  }

  double rainfall_eval(double u) const {
    return rainfall_spline.eval(u);
  }

  std::vector<double> rainfall_eval_range(std::vector<double> const& u) const {
    return rainfall_spline.r_eval(u);
  }

  /* if xi is not larger than every other x in the spline, recomputing the spline 
  will throw an error and will need to clear_points() and reconstruct spline */
  void rainfall_add_point(double xi, double yi) {
    rainfall_spline.add_point(xi, yi);
  }

  // safe add point
  void rainfall_add_point_sorted(double xi, double yi) {
    rainfall_spline.add_point_sorted(xi, yi);
  }

  /* not sure if these two are necessary for the rainfall model. could also move these
  to be general interpolator.h functions? */
  void rainfall_add_points(std::vector<double> const& x, std::vector<double> const& y) {
    util::check_length(x.size(), y.size());
    for (size_t i = 0; i < x.size(); ++i) {
      rainfall_add_point(x[i], y[i]); // refactor?
    }
  }

  // rainfall_add_points_sorted ?

  void rainfall_clear_points() {
    rainfall_spline.clear();
  }

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


  virtual void compute_rates() {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.set_rate(i, rainfall_eval(time) / (i+1));
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
  
private:
  interpolator::Interpolator rainfall_spline;
};

inline Rcpp::NumericMatrix get_state(const FF16_Environment environment) {
  return get_state(environment.canopy);
}


}

#endif
