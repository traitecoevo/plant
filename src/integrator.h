// -*-c++-*-
#ifndef TREE_INTEGRATOR_H_
#define TREE_INTEGRATOR_H_

#include <gsl/gsl_integration.h>
#include "functor.h"

namespace util {

// What is the pattern for use?  It seems that setting this up with
// some parameters (tolerances, bounds, etc) means we can then pass it
// around.  Use it as 
//    obj
class Integrator {
public:
  Integrator(double atol, double rtol, size_t max_iterations);
  ~Integrator();
  double integrate(DFunctor *f, double x_min, double x_max);
  
private:
  // Prevent copying by declaring these private with no definition.
  //Integrator(const Integrator &other);
  //Integrator& operator=(const Integrator &rhs);

  // TODO: atol/rtol are used throughout, but are abbreviations.
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
