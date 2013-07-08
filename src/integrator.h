// -*-c++-*-
#ifndef TREE_INTEGRATOR_H_
#define TREE_INTEGRATOR_H_

#include <gsl/gsl_integration.h>
#include "functor.h"

namespace util {

class Integrator {
public:
  Integrator(double atol_, double rtol_, size_t max_iterations_,
	     int quadrature_rule);
  ~Integrator();
  double integrate(DFunctor *f, double x_min, double x_max);
  static int gsl_rule(std::string rule);
  
private:
  // Prevent copying and assignment to prevent issues with gsl
  // pointers (see Spline for proper solution).
  Integrator(const Integrator &other);
  Integrator& operator=(Integrator other);

  double atol, rtol;
  size_t max_iterations;
  int quadrature_rule;
  bool with_singularities;

  double last_error;

  gsl_integration_workspace *workspace;
  gsl_function target_data;
};

namespace test {

double test_integrator(std::vector<double> pars, 
		       double x_min, double x_max,
		       std::string quadrature_rule);

}

}

#endif
