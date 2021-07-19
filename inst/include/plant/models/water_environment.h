// -*-c++-*-
#ifndef PLANT_PLANT_WATER_ENVIRONMENT_H_
#define PLANT_PLANT_WATER_ENVIRONMENT_H_

#include <plant/models/ff16_environment.h>

using namespace Rcpp;

namespace plant {

class Water_Environment : public FF16_Environment {
public:

  // Reuse methods from FF16_Environment
  using FF16_Environment::FF16_Environment;

  // Default constructor: three depths
  Water_Environment() {
    time = NA_REAL;
    canopy = Canopy();
    inflow_rate = 1.0;
    vars = Internals(3);
    set_soil_water_state(0.0);
  };

  // Custom constructor
  Water_Environment(Control control) {
    time = 0.0;
    canopy = Canopy(control);
    inflow_rate = control.soil_infiltration_rate;
    vars = Internals(control.soil_number_of_depths);
    set_soil_water_state(0.0);
  };

  double inflow_rate;

  // Patch interface
  virtual void compute_rates() {
    for (int i = 0; i < vars.state_size; i++) {
      vars.set_rate(i, inflow_rate / (i+1));
    }
  }

  // R interface
  void set_soil_water_state(double v) {
    for (int i = 0; i < vars.state_size; i++) {
      vars.set_state(i, v / (i+1));
    }
  }

  Internals r_internals() const { return vars; }
};

}

#endif
