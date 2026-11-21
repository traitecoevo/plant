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

// Writes `y` into row `i` of `out`, from column `offset`.
inline void set_row(Rcpp::NumericMatrix& out, size_t i,
                    const std::vector<double>& y, size_t offset = 0) {
  for (size_t j = 0; j < y.size(); ++j) {
    out(i, j + offset) = y[j];
  }
}

// The two ways run() grows a single plant from R. Each builds a Solver at
// odelia's default OdeControl() on a copy of `individual` in `environment`, runs
// it, and returns what run() reports; the plant passed in is left as it was.

// Grows the plant through `times`, which must be sorted and not before the
// environment's time, landing on each. Returns `state`, one row per time, and
// `individual`, the plant at each.
template <typename T, typename E>
Rcpp::List grow_to_times(const Individual<T,E>& individual,
                         const E& environment,
                         const std::vector<double>& times) {
  using System = IndividualRunner<T,E>;
  odelia::ode::Solver<System> solver(System(individual, environment),
                                     odelia::ode::OdeControl());
  solver.set_collect(false);
  const size_t n = times.size(), m = System::ode_size();
  Rcpp::NumericMatrix state(n, m);
  Rcpp::List at_time(n);
  for (size_t i = 0; i < n; ++i) {
    if (times[i] < solver.time()) {
      util::stop("times must be sorted and not before the environment's time, " +
                 util::to_string(solver.time()));
    }
    if (times[i] > solver.time()) {
      solver.advance_adaptive({solver.time(), times[i]});
    }
    set_row(state, i, solver.state());
    at_time[i] = Rcpp::wrap(solver.get_system_ref().individual);
  }
  return Rcpp::List::create(Rcpp::Named("state") = state,
                            Rcpp::Named("individual") = at_time);
}

// Grows the plant until ODE state `index` (from zero) reaches each of `sizes`
// (sorted) in turn, or to `time_max`. A size equal to the one before, or
// already reached, is reported where the run stands. Otherwise the time is
// located on the interpolant of the step that crosses the size, and the run
// lands there and goes on from there (odelia's Solver::advance_until), so a
// size's time, state and individual describe the same instant, integrated,
// not interpolated. Under Cash-Karp the interpolant is cubic Hermite, so the
// time is good to the interpolant's order, one below the step's: the state's
// entry `index` then differs from the size by the rate times that error, not
// by the location's tolerance.
//
// Returns `time` (NA for a size not reached), `state`, `individual` (NULL for
// a size not reached), `trajectory` (the start, then every step kept, as time
// and state), `timed_out` (time_max came before a size), `error` (the message
// where a step failed, else NULL) and `last`, the plant where the run stopped.
template <typename T, typename E>
Rcpp::List grow_to_sizes(const Individual<T,E>& individual,
                         const E& environment,
                         const std::vector<double>& sizes, size_t index,
                         double time_max) {
  using System = IndividualRunner<T,E>;
  const size_t n = sizes.size(), m = System::ode_size();
  if (index >= m) {
    util::stop("index " + util::to_string(index) +
               " is outside the ODE state, of size " + util::to_string(m));
  }
  odelia::ode::Solver<System> solver(System(individual, environment),
                                     odelia::ode::OdeControl());
  const double t_start = solver.time();
  const std::vector<double> y_start = solver.state();

  Rcpp::NumericVector time(n, NA_REAL);
  Rcpp::NumericMatrix state(n, m);
  std::fill(state.begin(), state.end(), NA_REAL);
  Rcpp::List at_size(n);
  bool timed_out = false;
  SEXP error = R_NilValue;

  // Locating on the interpolant evaluates no rates, so a tight tolerance
  // costs little; the interpolant's own error is the larger one.
  const odelia::root::tolerance tol{1e-10, 1e-10};
  for (size_t i = 0; i < n; ++i) {
    double t;
    if (i > 0 && sizes[i] == sizes[i - 1]) {
      t = time[i - 1];
    } else if (solver.state()[index] >= sizes[i]) {
      // Reached already: where a run starts on the size, or lands a hair past
      // one barely above the last.
      t = solver.time();
    } else {
      const double size = sizes[i];
      const auto reached = [index, size](double, const auto& y) {
        return y[index] - size;
      };
      std::optional<double> found;
      try {
        found = solver.advance_until(reached, time_max, tol);
      } catch (const std::exception& e) {
        error = Rcpp::wrap(std::string(e.what()));
        break;
      }
      if (!found) {
        timed_out = true;
        break;
      }
      t = *found;
    }
    time[i] = t;
    set_row(state, i, solver.state());
    at_size[i] = Rcpp::wrap(solver.get_system_ref().individual);
  }

  const size_t n_steps = solver.get_history_size();
  Rcpp::NumericMatrix trajectory(n_steps + 1, m + 1);
  trajectory(0, 0) = t_start;
  set_row(trajectory, 0, y_start, 1);
  std::vector<double> y(m);
  for (size_t k = 0; k < n_steps; ++k) {
    const System s = solver.get_history_step(k);
    s.ode_state(y.begin());
    trajectory(k + 1, 0) = s.ode_time();
    set_row(trajectory, k + 1, y, 1);
  }

  return Rcpp::List::create(
    Rcpp::Named("time") = time, Rcpp::Named("state") = state,
    Rcpp::Named("individual") = at_size,
    Rcpp::Named("trajectory") = trajectory,
    Rcpp::Named("timed_out") = timed_out, Rcpp::Named("error") = error,
    Rcpp::Named("last") = Rcpp::wrap(solver.get_system_ref().individual));
}

}
}

#endif /* PLANT_RUNNER */
