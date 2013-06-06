// -*-c++-*-
#ifndef TREE_LORENZ_H_
#define TREE_LORENZ_H_

#include <vector>
#include <Rcpp.h>

#include "ode_solver.h"

namespace ode {

namespace test {

class Lorenz {
public:
  Lorenz(double sigma, double R, double b);
  size_t size() const;
  void derivs(double time,
	      std::vector<double>::const_iterator y,
	      std::vector<double>::iterator dydt);

  void set_ode_state(std::vector<double> y, double t);
  std::vector<double> ode_state() const;
  double get_time() const;
  void step();
  void step_fixed(double step_size);
  void step_to(double time);
  void advance(double time_max);
  
  // * R interface
  std::vector<double> r_derivs(double time, 
			       std::vector<double> y);
  Rcpp::NumericMatrix r_run(std::vector<double> times,
			    std::vector<double> y);

private:
  const double sigma, R, b;
  Solver<Lorenz> solver;
  static const int ode_dimension = 3;
};

}

}

#endif
