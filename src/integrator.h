// -*-c++-*-
#ifndef TREE_INTEGRATOR_H_
#define TREE_INTEGRATOR_H_

#include <gsl/gsl_integration.h>
#include "functor.h"
#include <map>

namespace util {

class Integrator {
public:
  Integrator(double atol_, double rtol_, size_t max_iterations_,
	     int quadrature_rule);
  ~Integrator();
  Integrator(const Integrator &other);
  Integrator& operator=(Integrator other);

  double integrate(DFunctor *f, double x_min, double x_max);

  static int gsl_rule(std::string rule);
  static std::string gsl_rule_name(int rule);
  
private:

  typedef std::map<std::string, int> rules_type;
  static rules_type gsl_rule_table();

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
