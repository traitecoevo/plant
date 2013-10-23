// -*-c++-*-
#ifndef TREE_INTEGRATION_H_
#define TREE_INTEGRATION_H_

#include <Rcpp.h>
#include <list>

#include "functor.h"
#include "quadrature.h"
#include "integration_workspace.h"

namespace integration {

// This is the "QAG" algorithm from QUADPACK -- quadrature, adaptive,
// Gaussian.  Does not handle infinite intervals or singularities.
class QAG {
public:
  QAG(size_t rule, size_t max_iterations, double atol, double rtol);
  double integrate(util::DFunctor *f, double a, double b);

  // TODO: last_area inconsistent with QK::get_last_result()
  double get_last_area()       const {return area;}
  double get_last_error()      const {return error;}
  size_t get_last_iterations() const {return iteration;}

  // * R interface
  double r_integrate(util::RFunctionWrapper fun, double a, double b);

private:
  internal::workspace::point do_integrate(util::DFunctor *f,
					  double a, double b);
  bool initialise(util::DFunctor *f, double a, double b);
  bool refine(util::DFunctor *f);
  bool subinterval_too_small(double a1, double mid, double b2) const;

  QK q;
  internal::workspace w;

  // Control parameters
  size_t limit;
  double epsabs;
  double epsrel;

  // Intermediates
  double area, error;
  size_t iteration;
  size_t roundoff_type1, roundoff_type2;
};

}

#endif
