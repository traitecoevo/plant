// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_

#include <plant/control.h>
#include <plant/disturbance.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/ode_interface.h>
#include <plant/util.h>
#include <Rcpp.h>

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

  void compute_rates();

  // Dummy iterators: do nothing if the environment has no state.
  ode::iterator ode_state(ode::iterator it) const { return it; }
  ode::const_iterator set_ode_state(ode::const_iterator it) { return it; }
  ode::iterator ode_rates(ode::iterator it) const { return it; }
  int ode_size() const { return 0; }

  // Reset the environment.
  void clear() {
    time = 0.0;
    clear_environment();
  }

  void clear_environment() {}

  double offspring_arriving_dt() const {
    if (birth_rate.empty()) {
      Rcpp::stop("Cannot get offspring arrivals for empty environment");
    }
    return birth_rate[species_arriving_index];
  }

  void set_species_arriving_index(size_t x) {
    species_arriving_index = x;
  }

  // * R interface
  void r_set_species_arriving_index(util::index x) {
    set_species_arriving_index(x.check_bounds(birth_rate.size()));
  }

  void r_init_interpolators(const std::vector<double>& state) {}

  double get_environment_at_height(double height) { return 0.0; };

  double time;

  std::vector<double> birth_rate;
  size_t species_arriving_index;
};

}

#endif
