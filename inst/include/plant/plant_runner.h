// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_RUNNER_H_
#define PLANT_PLANT_PLANT_RUNNER_H_

#include <plant/plant.h>
#include <plant/environment.h>

namespace plant {
namespace tools {

template <typename T, typename E>
class PlantRunner {
public:
  typedef T strategy_type;
  typedef E environment_type;
  PlantRunner(Plant<T,E> plant_, environment_type environment_)
    : plant(plant_), environment(environment_) {
    plant.compute_rates(environment);
  }

  static size_t ode_size() {return Plant<T,E>::ode_size();}
  
  double ode_time() const {return environment.time;}
  ode::const_iterator set_ode_state(ode::const_iterator it, double time) {
    it = plant.set_ode_state(it);
    environment.time = time;
    plant.compute_rates(environment);
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    return plant.ode_state(it);
  }
  ode::iterator ode_rates(ode::iterator it) const {
    return plant.ode_rates(it);
  }
  
  Plant<T,E> plant;
  environment_type environment;
};

}
}

#endif
