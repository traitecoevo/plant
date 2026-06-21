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
    // Match FF16: loosen the light-availability spline tolerance from the
    // ResourceSpline default (1e-6) to 1e-4 for speed. The spline is rebuilt
    // every ODE step, so its construction dominates K93 runtime; 1e-6 was 100x
    // tighter than FF16 for no comparable accuracy need.
    light_availability = ResourceSpline(
        1e-4, // light_availability_spline_tol
        17,   // light_availability_spline_nbase
        16,   // light_availability_spline_max_depth
        true  // light_availability_spline_rescale_usually
    );
  };

  // Light interface
  ResourceSpline light_availability;

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

  virtual Rcpp::List r_get_state() const
  {
    return Rcpp::List::create(_["light_availability"] = time); //      light_availability);
  }

  // Core functions
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of light_availability spline
    // Note: extinction coefficient was already applied in strategy, so
    // f_compute_competition gives sum of projected leaf area (k L) across species. Just need to apply Beer's law, E = exp(- (k L))
    auto f_light_availability = [&](double height) -> double
    { return exp(-f_compute_competition(height)); };

    // Calculates the light_availability spline, by fitting to the function
    // `f_compute_competition` as a function of height

    light_availability.compute_environment(f_light_availability, height_max, rescale);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }

};

}

#endif
