// -*-c++-*-
#ifndef ADAPTIVESPLINE_H_
#define ADAPTIVESPLINE_H_

#include <list>
// #include <cmath>
#include <Rcpp.h>

#include "Spline.h"

class AdaptiveSpline : public Spline {
  // The target function takes two arguments: 'x' (the point to evaluate
  // a spline at) and '*data' a pointer to data used in the evaluation
  // of the target function.
  typedef double (*ASFun)(double x, void *data);
public:
  AdaptiveSpline();
  void set_bounds(double a_, double b_);
  void set_target(ASFun target_, void *data_);
  void set_control(double atol_, double rtol_, 
		   int nbase_, int max_depth_);
  bool construct_spline();

private:
  bool check_err(double y_true, double y_pred) const;
  void compute_spline();
  bool refine();
  void reset();

  void init();      // disable from base.
  // void init_self(); // disable from base?  Harmless enough.

  void check_bounds();

  double a, b, atol, rtol, dx, dxmin;
  int nbase, max_depth;
  ASFun target;
  void *data;

  // In contrast to the Spline's x and y, which are vectors, these are
  // lists so that we can easily add points in the middle of them.
  // 'yr' and 'yp' are y "real" and y "predicted", respectively.
  std::list<double> xx, yy;
  std::list<bool> zz;
};

// Using fixed method name 'target' so you can do
//   helper_spline<Foo>;
// to generate the correct function.
template <typename T>
double helper_spline(double x, void *data) {
  T *obj = static_cast<T*>(data);
  return (obj->target)(x);
}

// Using a method pointer as the second argument, so use is:
//   helper_spline<Foo, &Foo::bar>
// (assuming a method pointer of 'bar' here)
template <typename T, double (T::*Target)(double)>
double helper_spline(double x, void *data) {
  T *obj = static_cast<T*>(data);
  return (obj->*Target)(x);
}


#endif
