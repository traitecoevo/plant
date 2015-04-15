// -*-c++-*-
#ifndef TREE_FUNCTOR_H_
#define TREE_FUNCTOR_H_

// Something much flasher than this exists in Boost already, but we're
// reluctant to depend on boost because it's so massive.

#include <vector> // for test_functor
#include <Rcpp.h> // for SEXP

namespace util {

// Abstract class that defines a functor; the operator takes a single
// double as an argument and returns a single double as a result.
class DFunctor {
public:
  virtual double operator()(double) = 0;
  virtual ~DFunctor();
};

template <class T, double (T::*target)(double)> 
class Functor : public DFunctor {
public:
  Functor(T *obj_) : obj(obj_) {}
  virtual double operator()(double x) {
    return (obj->*target)(x);
  }
  
private:
  T* obj;
};

// This is a wrapper around the above for GSL.  It means you can set
// the "function" parameter as this function (util::helper_functor)
// and the "data" parameter as the pointer to the functor object.
double helper_functor(double x, void *data);

// This is a wrapper for creating a functor from an R function.
class RFunctionWrapper : public DFunctor {
public:
  RFunctionWrapper(Rcpp::Function fun_);
  double target(double x);
  double operator()(double x);
private:
  Rcpp::Function fun;
};

// To prepare for the integration in `compute_assimilation` we need to
// convert the function `compute_assimilation_x(double, util::Spline*)
// to take just a double as an argument.  Boost has the ability to
// bind arguments which would be nice here, but we're avoiding
// depending on that for the time being.
//
// This binds the second argument, assuming a const method.  All a bit
// of a hack, but it does seem to work correctly.
template <class T, class T2, double (T::*target)(double, T2) const>
class FunctorBind1 : public util::DFunctor {
public:
  FunctorBind1(const T *obj_, T2 arg2_) : obj(obj_), arg2(arg2_) {}
  virtual double operator()(double x) {
    return (obj->*target)(x, arg2);
  }

private:
  const T* obj;
  T2 arg2;
};

// Painfully copied from Plant, but because we need a non-const
// version.
template <class T, class T2, double (T::*target)(double, T2)>
class FunctorBind2 : public util::DFunctor {
public:
  FunctorBind2(T *obj_, T2 arg2_) : obj(obj_), arg2(arg2_) {}
  virtual double operator()(double x) {
    return (obj->*target)(x, arg2);
  }
private:
  T* obj;
  T2 arg2;
};

// This is used for testing in a couple of places.
namespace test {

class Quadratic {
public:
  Quadratic(double a_, double b_, double c_): a(a_), b(b_), c(c_) {}
  // Can be const, except that find_root.cpp needs it not to be.
  double mytarget(double x) { return (a * x + b) * x + c; }
private:
  const double a, b, c;
};

std::vector<double> test_functor(std::vector<double> x, 
				 std::vector<double> pars);
}

}

RCPP_EXPOSED_CLASS_NODECL(util::RFunctionWrapper)

#endif
