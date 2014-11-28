#include <tree/find_root.h>

#include <tree/util.h>

namespace util {

RootFinder::RootFinder(double atol_, double rtol_, int max_iterations_)
  : atol(atol_),
    rtol(rtol_),
    max_iterations(max_iterations_),
    iterations(0),
    last_error(NA_REAL) {
  F.function = &helper_functor;
  solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
}

RootFinder::~RootFinder() {
  gsl_root_fsolver_free(solver);
}

double RootFinder::root(DFunctor *f, 
			double x_min, double x_max) {
  F.params = static_cast<void*>(f);
  gsl_root_fsolver_set(solver, &F, x_min, x_max);
  iterations = 0;

  int status = GSL_CONTINUE;
  double r;
  do {
    iterations++;
    gsl_root_fsolver_iterate(solver);
    r = gsl_root_fsolver_root(solver);
    x_min = gsl_root_fsolver_x_lower(solver);
    x_max = gsl_root_fsolver_x_upper(solver);
    status = gsl_root_test_interval(x_min, x_max, atol, rtol);
  } while (status == GSL_CONTINUE && iterations < max_iterations);
  last_error = x_max - x_min;

  if (status != GSL_SUCCESS)
    Rcpp::stop("RootFinder failure: status " +
	       util::to_string(status) + " in " +
	       util::to_string(iterations) + " (err: " +
	       util::to_string(last_error) + ")");

  return r;
}

namespace test {

double test_find_root(std::vector<double> pars,
		      double x_min, double x_max) {
  util::check_length(pars.size(), 3);
  Quadratic obj(pars[0], pars[1], pars[2]);
  Functor<test::Quadratic, &test::Quadratic::mytarget> fun(&obj);

  // Control parameters for the root solver
  const double atol = 1e-6, rtol = 1e-6;
  const int max_iterations = 1000;
  RootFinder root(atol, rtol, max_iterations);
  return root.root(&fun, x_min, x_max);
}

double test_find_value(std::vector<double> pars, double value,
		       double x_min, double x_max) {
  util::check_length(pars.size(), 3);
  Quadratic obj(pars[0], pars[1], pars[2]);
  FunctorRoot<test::Quadratic, &test::Quadratic::mytarget> 
    fun(&obj, value);

  // Control parameters for the root solver
  const double atol = 1e-6, rtol = 1e-6;
  const int max_iterations = 1000;
  RootFinder root(atol, rtol, max_iterations);
  return root.root(&fun, x_min, x_max);
}

}

}
