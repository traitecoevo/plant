// -*-c++-*-
#ifndef TREE_FIND_ROOT_
#define TREE_FIND_ROOT_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "functor.h"

namespace util {

double find_root(DFunctor *f, double x_lo, double x_hi);

namespace test {
double test_find_root(std::vector<double> pars, double x_lo, double x_hi);
}

}

#endif
