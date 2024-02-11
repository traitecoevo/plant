// Built from  inst/include/plant/models/ff16_environment.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_ENVIRONMENT_H_
#define PLANT_PLANT_TF24_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <plant/interpolator.h>

using namespace Rcpp;

namespace plant {

class TF24_Environment : public Environment {
public:
  // constructor for R interface - default settings can be modified
  // except for soil_number_of_depths and light_availability_spline_rescale_usually
  // which are only updated on construction
  TF24_Environment(bool light_availability_spline_rescale_usually = false,
                   int soil_number_of_depths = 0) {
    time = 0.0;
    
    light_availability = ResourceSpline();
    light_availability.spline_rescale_usually = light_availability_spline_rescale_usually;

    vars = Internals(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, 0.0));
  };

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline light_availability;

  // Ability to prescribe a fixed value
  // TODO: add setting to set other variables like water
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

  // Pre-compute resources available in the environment, as a function of height
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


inline Rcpp::List get_state(const TF24_Environment environment, double time) {
  auto ret = get_state(environment.extrinsic_drivers, time);
  ret["light_availability"] = get_state(environment.light_availability);
  return ret;
}
}

#endif
