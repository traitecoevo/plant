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
  using value_type = typename T::value_type;

  IndividualRunner(Individual<T,E> individual_, E environment_)
    : individual(individual_), environment(environment_) {
    individual.compute_rates(environment);
  }

  static size_t ode_size() {return Individual<T,E>::ode_size();}
  
  double ode_time() const {return environment.time;}
  template <typename It> It set_ode_state(It it, double time) {
    it = individual.set_ode_state(it);
    environment.time = time;
    individual.compute_rates(environment);
    return it;
  }
  template <typename It> It ode_state(It it) const {
    return individual.ode_state(it);
  }
  template <typename It> It ode_rates(It it) const {
    return individual.ode_rates(it);
  }
  
  Individual<T,E> individual;
  E environment;
};

}
}

#endif /* PLANT_RUNNER */
