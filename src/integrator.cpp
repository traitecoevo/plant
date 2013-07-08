#include "integrator.h"

#include "util.h"

namespace util {

Integrator::Integrator(double atol_, double rtol_,
		       size_t max_iterations_, int quadrature_rule_)
  : atol(atol_),
    rtol(rtol_),
    max_iterations(max_iterations_),
    quadrature_rule(quadrature_rule_),
    with_singularities(quadrature_rule == TREE_GSL_INTEG_QAGS),
    nonadaptive(quadrature_rule == TREE_GSL_INTEG_NONADAPTIVE),
    last_error(NA_REAL) {
  target_data.function = &helper_functor;
  workspace = gsl_integration_workspace_alloc(max_iterations);
}

Integrator::~Integrator() {
  gsl_integration_workspace_free(workspace);
}

Integrator::Integrator(const Integrator& other)
  : atol(other.atol),
    rtol(other.rtol),
    max_iterations(other.max_iterations),
    quadrature_rule(other.quadrature_rule),
    with_singularities(other.with_singularities),
    nonadaptive(other.nonadaptive),
    last_error(other.last_error) {
  target_data.function = &helper_functor;
  workspace = gsl_integration_workspace_alloc(max_iterations);
}

Integrator& Integrator::operator=(Integrator other) {
  using std::swap;
  swap(atol,               other.atol);
  swap(rtol,               other.rtol);
  swap(max_iterations,     other.max_iterations);
  swap(quadrature_rule,    other.quadrature_rule);
  swap(with_singularities, other.with_singularities);
  swap(nonadaptive,        other.nonadaptive);
  swap(last_error,         other.last_error);
  swap(workspace,          other.workspace);
  swap(target_data,        other.target_data);
  return *this;
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
  } else if (nonadaptive) {
    size_t n_eval;
    gsl_integration_qng(&target_data, x_min, x_max,
			atol, rtol,
			&result, &last_error, &n_eval);
  } else {
    gsl_integration_qag(&target_data, x_min, x_max,
			atol, rtol, max_iterations, quadrature_rule,
			workspace,
			&result, &last_error);
  }
  return result;
}

int Integrator::gsl_rule(std::string rule) {
  const rules_type rules = gsl_rule_table();
  rules_type::const_iterator it = rules.find(rule);
  if (it == rules.end())
    ::Rf_error("Unknown rule");
  return it->second;
}

std::string Integrator::gsl_rule_name(int rule) {
  const rules_type rules = gsl_rule_table();
  rules_type::const_iterator it = rules.begin();
  while (it != rules.end()) {
    if (it->second == rule)
      return it->first;
    ++it;
  }
  ::Rf_error("Unknown rule");
}

// I don't see why this is necessary -- it only matters for the method
// declaration, and the declaration for 'rules' and within the class
// definition pass without error.
typedef std::map<std::string, int> rules_type;
rules_type Integrator::gsl_rule_table() {
  rules_type rules;
  rules["GAUSS15"] = GSL_INTEG_GAUSS15;
  rules["GAUSS21"] = GSL_INTEG_GAUSS21;
  rules["GAUSS31"] = GSL_INTEG_GAUSS31;
  rules["GAUSS41"] = GSL_INTEG_GAUSS41;
  rules["GAUSS51"] = GSL_INTEG_GAUSS51;
  rules["GAUSS61"] = GSL_INTEG_GAUSS61;
  rules["QAGS"]    = -1;
  rules["NONADAPTIVE"] = -2;
  return rules;
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
