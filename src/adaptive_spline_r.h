// -*-c++-*-
#ifndef TREE_ADAPTIVE_SPLINE_R_H_
#define TREE_ADAPTIVE_SPLINE_R_H_

#include <Rcpp.h>

#include "adaptive_spline.h"
#include "functor.h"

namespace spline {

namespace test {

Spline test_adaptive_spline(SEXP fun, SEXP env,
			    double a, double b,
			    bool akima_spline);
}


}

#endif
