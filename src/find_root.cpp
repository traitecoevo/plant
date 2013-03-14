#include <Rcpp.h>

#include "find_root.h"

namespace util {

double find_root(DFunctor *f, double x_lo, double x_hi) {
  double atol = 0.0, rtol = 0.001;
  int iter = 0, max_iter = 100;

  gsl_function F;
  F.function = &helper_functor;
  F.params = (void*)(f);
  
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  int status;
  double r;
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, atol, rtol);
  } while (status == GSL_CONTINUE && iter < max_iter);
  
  // clean up
  gsl_root_fsolver_free(s);

  // TODO: Detect false convergence and throw.
  // Rprintf("Converged with status %d in %d iterations (err: %2.5f)\n", 
  // 	  status, iter, x_hi - x_lo);


  return r;
}

namespace test {

double test_find_root(std::vector<double> pars, 
		      double x_lo, double x_hi) {
  if ( (int)pars.size() != 3 )
    Rf_error("Expected parameters of length 3");
  Quadratic obj(pars[0], pars[1], pars[2]);
  
  Functor<test::Quadratic, &test::Quadratic::mytarget> fun(&obj);
  return find_root(&fun, x_lo, x_hi);
}

}

}
