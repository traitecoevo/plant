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

// TODO: setting `key` to different values gives different integration
// rules:
//   GSL_INTEG_GAUSS15
//   GSL_INTEG_GAUSS21
//   GSL_INTEG_GAUSS31
//   GSL_INTEG_GAUSS41
//   GSL_INTEG_GAUSS51
//   GSL_INTEG_GAUSS61
// corresponding to the 15, 21, 31, 41, 51 and 61 point
// Gauss-Kronrod rules.  The higher-order rules give better accuracy
// for smooth functions, while lower-order rules save time when the
// function contains local difficulties, such as discontinuities.
// 
// Might be worth tweaking this.  It can probably be sorted out at the
// construction phase with the other control parameters.  It might be
// possible to do some `Lookup`-like control for use from within R,
// but only if it is also straightforward to change from C++.
double Integrator::integrate(DFunctor *f, double x_min, double x_max) {
  target_data.params = f;
  double result;

  // Not sure if we have a singularity, but we'll want to swich
  // between these two algorithms if we do.
  // gsl_integration_qag(&target_data, x_min, x_max, 
  // 		      atol, rtol, max_iterations,
  // 		      GSL_INTEG_GAUSS15, workspace, 
  // 		      &result, &last_error);
  // In Rmath, we were using Rdqags, which should be the same
  // algorithm as qags.
  // 
  // However, we generally know the points in advance for the problems
  // that I've been having -- they are one of the end points.

  gsl_integration_qags(&target_data, x_min, x_max, 
		       atol, rtol, max_iterations,
		       workspace, 
		       &result, &last_error);
  // We do not do anything with the return status, as GSL will already
  // throw an error for us if the integration fails.  It does not
  // appear to be possible to report the number of iterations used.
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
