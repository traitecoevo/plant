#include "integrator.h"

#include "util.h"

namespace util {

Integrator::Integrator(double atol, double rtol, size_t max_iterations) 
  : atol(atol), rtol(rtol), max_iterations(max_iterations) {
  target_data.function = &helper_functor;
  workspace = gsl_integration_workspace_alloc(max_iterations);
}

Integrator::~Integrator() {
  gsl_integration_workspace_free(workspace);
}

// NOTE: We do not do anything with the return status, as GSL will
// already throw an error for us if the integration fails.  It does
// not appear to be possible to report the number of iterations used.
double Integrator::integrate(DFunctor *f, double x_min, double x_max) {
  target_data.params = f;
  double result = 0.0;
  gsl_integration_qags(&target_data, x_min, x_max, 
		       atol, rtol, max_iterations,
		       workspace, 
		       &result, &last_error);
  return result;
}

namespace test {

double test_integrator(std::vector<double> pars, 
		       double x_min, double x_max) {
  // Set up the usual quadratic problem:
  util::check_length(pars.size(), 3);
  Quadratic obj(pars[0], pars[1], pars[2]);
  Functor<test::Quadratic, &test::Quadratic::mytarget> fun(&obj);

  // Control parameters for the integrator
  const double atol = 1e-6, rtol = 1e-6;
  const size_t max_iterations = 1000;

  Integrator integrator(atol, rtol, max_iterations);
  return integrator.integrate(&fun, x_min, x_max);
}

}

}
