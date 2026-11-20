// -*-c++-*-
#ifndef PLANT_RUNNER
#define PLANT_RUNNER

#include <plant/community/individual.h>
#include <plant/models/environment.h>
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

// Steps `solver` until ODE state `index` (from zero) reaches `size`, or to
// `time_max`, and returns the time it did, or NA where time_max came first.
// The time is located on the interpolant of the step that crossed `size`, and
// the solver then steps from that step's start to land on it (odelia's
// Solver::advance_until), so the solver's state and the individual it holds
// are the plant at the returned time, integrated, not interpolated. Under
// Cash-Karp the interpolant is cubic Hermite, so the time is good to the
// interpolant's order, one below the step's: the state's entry `index` then
// differs from `size` by the rate times that error, not by the location's
// tolerance.
template <typename Solver>
double advance_to_size(Solver& solver, size_t index, double size,
                       double time_max) {
  using System = typename Solver::system_type;
  if (index >= System::ode_size()) {
    util::stop("index " + util::to_string(index) +
               " is outside the ODE state, of size " +
               util::to_string(System::ode_size()));
  }
  const auto reached = [index, size](double, const auto& y) {
    return y[index] - size;
  };
  // Locating on the interpolant evaluates no rates, so a tight tolerance
  // costs little; the interpolant's own error is the larger one.
  const odelia::root::tolerance tol{1e-10, 1e-10};
  const std::optional<double> t = solver.advance_until(reached, time_max, tol);
  return t ? *t : NA_REAL;
}

// The steps `solver` has kept, one row each: the time, then the ODE state.
// These are the steps that stand, so a step that crossed a size and was taken
// back by advance_to_size() is not among them; the step that landed on the
// size is.
template <typename Solver>
Rcpp::NumericMatrix steps_taken(const Solver& solver) {
  using System = typename Solver::system_type;
  const size_t n = solver.get_history_size();
  const size_t m = System::ode_size();
  Rcpp::NumericMatrix out(n, m + 1);
  std::vector<double> y(m);
  for (size_t i = 0; i < n; ++i) {
    const System s = solver.get_history_step(i);
    s.ode_state(y.begin());
    out(i, 0) = s.ode_time();
    for (size_t j = 0; j < m; ++j) {
      out(i, j + 1) = y[j];
    }
  }
  return out;
}

}
}

#endif /* PLANT_RUNNER */
