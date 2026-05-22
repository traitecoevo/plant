// -*-c++-*-
#ifndef PLANT_RUNNER
#define PLANT_RUNNER

#include <plant/individual.h>
#include <plant/environment.h>

namespace plant {
namespace tools {

template <typename T, typename E>
class IndividualRunner {
public:
  IndividualRunner(Individual<T,E> individual_, E environment_)
    : individual(individual_), environment(environment_) {
    individual.compute_rates(environment);
  }

  static size_t ode_size() {return Individual<T,E>::ode_size();}
  
  double ode_time() const {return environment.time;}
  ode::const_iterator set_ode_state(ode::const_iterator it, double time) {
    it = individual.set_ode_state(it);
    environment.time = time;
    individual.compute_rates(environment);
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    return individual.ode_state(it);
  }
  ode::iterator ode_rates(ode::iterator it) const {
    return individual.ode_rates(it);
  }
  
  Individual<T,E> individual;
  E environment;
};

}
}

#endif /* PLANT_RUNNER */
