#include "integrator.h"

#include "util.h"

namespace util {

Integrator::Integrator(double atol_, double rtol_,
		       size_t max_iterations_, int quadrature_rule_)
  : atol(atol_),
    rtol(rtol_),
    max_iterations(max_iterations_),
    quadrature_rule(quadrature_rule_),
    with_singularities(quadrature_rule < 0),
    last_error(NA_REAL) {
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
  if (with_singularities) {
    gsl_integration_qags(&target_data, x_min, x_max,
			 atol, rtol, max_iterations,
			 workspace,
			 &result, &last_error);
  } else {
    gsl_integration_qag(&target_data, x_min, x_max,
			atol, rtol, max_iterations, quadrature_rule,
			workspace,
			&result, &last_error);
  }
  return result;
}

// NOTE: Static function:
int Integrator::gsl_rule(std::string rule) {
  int ret = 0;
  if      (rule == "GAUSS15") ret = GSL_INTEG_GAUSS15;
  else if (rule == "GAUSS21") ret = GSL_INTEG_GAUSS21;
  else if (rule == "GAUSS31") ret = GSL_INTEG_GAUSS31;
  else if (rule == "GAUSS41") ret = GSL_INTEG_GAUSS41;
  else if (rule == "GAUSS51") ret = GSL_INTEG_GAUSS51;
  else if (rule == "GAUSS61") ret = GSL_INTEG_GAUSS61;
  else if (rule == "QAGS")    ret = -1; // NOTE: hack
  else ::Rf_error("Unknown rule");
  return ret;
}

namespace test {

double test_integrator(std::vector<double> pars, 
		       double x_min, double x_max,
		       std::string quadrature_rule) {
  // Set up the usual quadratic problem:
  util::check_length(pars.size(), 3);
  Quadratic obj(pars[0], pars[1], pars[2]);
  Functor<test::Quadratic, &test::Quadratic::mytarget> fun(&obj);

  // Control parameters for the integrator
  const double atol = 1e-6, rtol = 1e-6;
  const size_t max_iterations = 1000;
  const int rule = Integrator::gsl_rule(quadrature_rule);

  Integrator integrator(atol, rtol, max_iterations, rule);
  return integrator.integrate(&fun, x_min, x_max);
}

}

}
