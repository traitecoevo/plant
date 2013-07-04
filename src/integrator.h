// -*-c++-*-
#ifndef TREE_INTEGRATOR_H_
#define TREE_INTEGRATOR_H_

#include <gsl/gsl_integration.h>
#include "functor.h"

namespace util {

class Integrator {
public:
  Integrator(double atol_, double rtol_, size_t max_iterations_);
  ~Integrator();
  double integrate(DFunctor *f, double x_min, double x_max);
  
private:
  // Prevent copying and assignment to prevent issues with gsl
  // pointers (see Spline for proper solution).
  Integrator(const Integrator &other);
  Integrator& operator=(Integrator other);

  double atol, rtol;
  size_t max_iterations;

  double last_error;

  gsl_integration_workspace *workspace;
  gsl_function target_data;
};

namespace test {

double test_integrator(std::vector<double> pars, 
		       double x_min, double x_max);

}

}

#endif
