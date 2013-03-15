// -*-c++-*-
#ifndef TREE_FIND_ROOT_
#define TREE_FIND_ROOT_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "functor.h"

namespace util {

// TODO: not totally clear what the best name is here for the class.
// Should be a noun (as is Spline, Solver, etc).
class RootFinder {
public:
  RootFinder(double atol, double rtol, int max_iterations);
  ~RootFinder();
  double root(DFunctor *f, double x_min, double x_max);

private:
  // Prevent copying and assignment
  // RootFinder(const RootFinder other&);
  // RootFinder& operator=(const RootFinder rhs=);

  // TODO: Here, and elsewhere, these are abbreviations.
  double atol, rtol;
  int max_iterations;
  int iterations;

  gsl_function F;
  gsl_root_fsolver *solver;

  // I had something like this last time.
  double last_error;
};

// This may be useful as more often than not we are *actually*
// interested in where the target function takes a particular value
// `value`, rather than when it is zero.
// 
// With this we can do
// FunctorRoot<MyClass,&MyClass::target> f(this, 2.3);
// and then pass &f as the functor to the root finder:
// 
// RootFinder root_finder(atol, rtol, max_iterations);
// root_finder.find_root(&f, x_min, x_max)
// 
// This could probaby be easier if we made a functioniod so that
// something like
// 
// root_finder.find_root(FunctorRoot<MyClass,&MyClass::target>(this, 2.3),
//                       x_min, x_max)
//
// worked, but that's actually not much clearer.

template <class T, double(T::*target)(double)>
class FunctorRoot : public DFunctor {
public:
  FunctorRoot(T *obj, double value) : obj(obj), value(value) {}
  virtual double operator()(double x) {
    return (obj->*target)(x) - value;
  }
private:
  T* obj;
  double value;
};
  
namespace test {
double test_find_root(std::vector<double> pars, double x_min, double x_max);
double test_find_value(std::vector<double> pars, double value,
		       double x_min, double x_max);
}

}

#endif
