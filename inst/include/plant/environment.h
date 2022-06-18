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
  ExtrinsicDrivers extrinsic_drivers;
};
}
#endif
