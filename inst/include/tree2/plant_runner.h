// -*-c++-*-
#ifndef TREE_PLANT_RUNNER_H_
#define TREE_PLANT_RUNNER_H_

#include <tree2/plant.h>
#include <tree2/environment.h>

namespace tree2 {
namespace tools {

struct PlantRunner {
  PlantRunner(Plant plant_, Environment environment_)
    : plant(plant_), environment(environment_) {
    plant.compute_vars_phys(environment);
  }
  static size_t ode_size() {return Plant::ode_size();}
  ode::const_iterator set_ode_state(ode::const_iterator it) {
    it = plant.set_ode_state(it);
    plant.compute_vars_phys(environment);
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    return plant.ode_state(it);
  }
  ode::iterator ode_rates(ode::iterator it) const {
    return plant.ode_rates(it);
  }
  Plant plant;
  Environment environment;
};

}
}

#endif
