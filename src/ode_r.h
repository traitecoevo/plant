// -*-c++-*-
#ifndef TREE_ODE_R_H_
#define TREE_ODE_R_H_

#include <Rcpp.h>
#include <vector>
#include <cstddef>

#include "ode_solver.h"

namespace ode {

class OdeR {
public:
  OdeR(SEXP fun_, SEXP env_, SEXP pars_);
  OdeR(SEXP fun_, SEXP env_, SEXP pars_, OdeControl control);
  size_t size() const;
  void derivs(double time, iterator_const y, iterator dydt);

  void set_ode_state(std::vector<double> y, double time);
  std::vector<double> ode_state() const;
  double get_time() const;
  std::vector<double> get_times() const;

  void step();
  void step_fixed(double step_size);
  void step_to(double time_max);
  void advance(double time_max);
  void advance_fixed(std::vector<double> times);

  void reset();

  // * R interface
  std::vector<double> r_derivs(double time, std::vector<double> y);
  Rcpp::NumericMatrix r_run(std::vector<double> times,
			    std::vector<double> y);
  OdeControl r_control() const;
  SEXP r_pars() const;

private:
  SEXP target(double time, SEXP y);

  SEXP fun, env, pars;
  Solver<OdeR> solver;
};

}

#endif
