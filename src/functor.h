// -*-c++-*-
#ifndef TREE_FUNCTOR_H_
#define TREE_FUNCTOR_H_

// Something much flasher than this exists in Boost already, but we're
// reluctant to depend on boost because it's so massive.

#include <vector> // for test_functor

namespace util {

// Abstract class that defines a functor; the operator takes a single
// double as an argument and returns a single double as a result.
// 
// TODO: Possibly this would be nicer to template as 
// template <class R, class A>
// class Functor {
// public:
//   virtual R operator()(A) = 0;
// };
class DFunctor {
public:
  virtual double operator()(double) = 0;
};

template <class T, double (T::*target)(double)> 
class Functor : public DFunctor {
public:
  Functor(T *obj) : obj(obj) {}
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

// This is used for testing in a couple of places.
namespace test {

class Quadratic {
public:
  Quadratic(double a, double b, double c): a(a), b(b), c(c) {}
  double mytarget(double x) { return (a * x + b) * x + c; }
private:
  const double a, b, c;
};

std::vector<double> test_functor(std::vector<double> x, 
				 std::vector<double> pars);
}

}

#endif
