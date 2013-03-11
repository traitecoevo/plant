// -*-c++-*-
#ifndef TREE_LORENZ_H_
#define TREE_LORENZ_H_

#include <vector>
#include <Rcpp.h>

#include "ode_solver.h"

namespace ode {

class Lorenz {
public:
  Lorenz(double sigma, double R, double b);
  void derivs(double time,
	      std::vector<double>::const_iterator y,
	      std::vector<double>::iterator dydt);
  unsigned int size() const { return 3; }
  
  std::vector<double> r_derivs(double time, 
			       std::vector<double> y);

  // These are going to be tedious to expose, but in this case we're
  // really not trying to do any abstraction at all.
  void ode_set_state(std::vector<double> y, double t) {
    solver.set_state(y, t);
  }
  std::vector<double> ode_get_state() const {
    return solver.get_state();
  }
  double ode_get_time() const {
    return solver.get_time();
  }
  void ode_step() {
    solver.step();
  }
  void ode_step_fixed(double step_size) { 
    solver.step_fixed(step_size); 
  }
  void ode_advance(double time_max) {
    solver.advance(time_max);
  }
  Rcpp::NumericMatrix ode_r_run(std::vector<double> times, 
				std::vector<double> y) {
    return solver.r_run(times, y);
  }

private:
  const double sigma, R, b;
  Solver<Lorenz> solver;
};

}

#endif
