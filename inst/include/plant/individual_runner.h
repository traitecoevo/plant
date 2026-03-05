// -*-c++-*-
#ifndef PLANT_RUNNER
#define PLANT_RUNNER

#include <plant/individual.h>
#include <plant/environment.h>
#include <odelia/ode_solver.hpp>

namespace plant {
namespace tools {

template <typename T, typename E>
class IndividualRunner {
public:
  using value_type = double;

  IndividualRunner(Individual<T,E> individual_, E environment_)
    : individual(individual_), environment(environment_) {
    individual.compute_rates(environment);
  }

  static size_t ode_size() {return Individual<T,E>::ode_size();}
  
  double ode_time() const {return environment.time;}
  odelia::ode::const_iterator set_ode_state(odelia::ode::const_iterator it, double time) {
    it = individual.set_ode_state(it);
    environment.time = time;
    individual.compute_rates(environment);
    return it;
  }
  odelia::ode::iterator ode_state(odelia::ode::iterator it) const {
    return individual.ode_state(it);
  }
  odelia::ode::iterator ode_rates(odelia::ode::iterator it) const {
    return individual.ode_rates(it);
  }
  
  Individual<T,E> individual;
  E environment;
};

}
}

#endif /* PLANT_RUNNER */
