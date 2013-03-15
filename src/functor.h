// -*-c++-*-
#ifndef TREE_FUNCTOR_
#define TREE_FUNCTOR_

// Something much flasher than this exists in Boost already, but we're
// reluctant to depend on boost because it's so massive.

#include <vector> // for test_functor

namespace util {

// Abstract class that defines a functor; the operator takes a single
// double as an argument and returns a single double as a result.
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

// But is that any easier than the approach with AdaptiveSpline of
// having
// template <typename T, double (T::*Target)(double)>
// double helper_function(double x, void *data) {
//   T *obj = static_cast<T*>(data);
//   return (obj->*Target)(x);
// }
// ?
// 
// It makes the interface very slightly simpler, so possibly.

// Bind a second argument to return a functor for the first.
// This seems a bit crap in terms of generalising and naming.  Boost
// would certainly be easier here.
template <class T, class T2, double (T::*target)(double, T2)>
class FunctorBind1 : public DFunctor {
public:
  FunctorBind1(T *obj, T2 arg2) : obj(obj), arg2(arg2) {}
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
