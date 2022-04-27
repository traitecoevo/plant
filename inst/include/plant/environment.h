// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_

#include <plant/control.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/ode_interface.h>
#include <plant/internals.h>
#include <plant/util.h>
#include <unordered_map>
#include <Rcpp.h>

#include <plant/extrinsic_drivers.h>

using namespace Rcpp;

namespace plant {

class Environment {
public:
  template <typename Function>
  void compute_environment(Function f, double height_max);
  template <typename Function>
  void rescale_environment(Function f, double height_max);

  void set_fixed_environment(double value, double height_max);
  void set_fixed_environment(double value);
  void init_interpolators(const std::vector<double>& state);

  // ODE interface: do nothing if the environment has no state.
  size_t ode_size() const { return vars.state_size; }
  virtual void compute_rates(std::vector<double> const& resource_depletion){};

  ode::const_iterator set_ode_state(ode::const_iterator it) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.states[i] = *it++;
    }
    return it;
  }

  ode::iterator ode_state(ode::iterator it) const {
    for (size_t i = 0; i < vars.state_size; i++) {
      *it++ = vars.states[i];
    }
    return it;
  }

  ode::iterator ode_rates(ode::iterator it) const {
    for (size_t i = 0; i < vars.state_size; i++) {
      *it++ = vars.rates[i];
    }
    return it;
  }

  // Reset the environment.
  void clear() {
    time = 0.0;
    clear_environment();
  }

  void clear_environment() {}

  void r_init_interpolators(const std::vector<double>& state) {}

  double get_environment_at_height(double height) { return 0.0; };
  
  virtual ~Environment() = default;

  double time;

  size_t species_arriving_index;

  Internals vars;

  /* EXSTRINSIC DRIVERS API*/

  // initialise spline of driver with x, y control points
  void set_extrinsic_driver(std::string driver_name, std::vector<double> const& x, std::vector<double> const& y) {
    // if we wanted to be faster we could skip this check (but less safe)
    if (extrinsic_drivers.count(driver_name) == 0) {
      util::stop(driver_name + " doesn't exist in the list of extrinsic_drivers.");
    } else {
      extrinsic_drivers.at(driver_name).init(x, y);
      extrinsic_driver_extrapolation(driver_name, false); // default no extrapolation
    }
  }

  void extrinsic_driver_extrapolation(std::string driver_name, bool extrapolate) {
    extrinsic_drivers.at(driver_name).set_extrapolate(extrapolate);
  }

  // evaluate/query interpolated spline for driver at point u, return s(u), where s is interpolated function
  double extrinsic_driver_evaluate(std::string driver_name, double u) const {
    return extrinsic_drivers.at(driver_name).eval(u);
  }

  // evaluate/query interpolated spline for driver at vector of points, return vector of values
  std::vector<double> extrinsic_driver_evaluate_range(std::string driver_name, std::vector<double> u) const {
    return extrinsic_drivers.at(driver_name).r_eval(u);
  }

  // returns the name of each active driver - useful for R output
  std::vector<std::string> get_extrinsic_driver_names() {
    auto ret = std::vector<std::string>();
    for (auto const& driver : extrinsic_drivers) {
      ret.push_back(driver.first);
    }
    return ret;
  }

protected:
  std::unordered_map<std::string, interpolator::Interpolator> extrinsic_drivers;
};

}

#endif
