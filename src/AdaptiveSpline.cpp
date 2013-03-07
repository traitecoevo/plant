#include <Rcpp.h>
#include "AdaptiveSpline.h"

#include <gsl/gsl_nan.h>

#include "util.h"


// Constructor that sets things to sensible defaults.  The bounds are
// impossible, and the target function is NULL, both of which will
// prevent construction of the spline without calling `set_bounds()`
// and `set_targets()`.  Setting the control parameters (via
// `set_control`) is optional, and the bounds should be reasonable for
// many uses.
AdaptiveSpline::AdaptiveSpline() : Spline() { 
  a = GSL_POSINF;
  b = GSL_NEGINF;

  target = NULL;

  atol = 1e-6;
  rtol = 1e-6;
  nbase = 17;
  max_depth = 8;
}

// TODO: Some translation here with R and GSL infinities might be
// useful.
void AdaptiveSpline::set_bounds(double a_, double b_) {
  a = a_;
  b = b_;
  check_bounds();
}

// TODO: This requires some effort to explain -- what is the target
// function and what is the data pointer?
void AdaptiveSpline::set_target(ASFun target_, void *data_) {
  target = target_;
  data = data_;
}

// Set *all* the control parameters.
// 
// TODO: This could be refined, as currently this requires that all
// are known to be set well.  It's not like there is any checking on
// these (positivity, etc).  Could make nbase and max_depth be
// unsigned integers?  Or could write a series of get/set.
void AdaptiveSpline::set_control(double atol_, double rtol_, 
				 int nbase_, int max_depth_) {
  atol = atol_;
  rtol = rtol_;
  nbase = nbase_;
  max_depth = max_depth_;
}

// This goes through and adaptively builds the spline.
bool AdaptiveSpline::construct_spline() {
  if ( target == NULL )
    Rf_error("Target function not set");
  check_bounds();
  reset();

  dx = (b - a) / (nbase - 1);
  dxmin = dx / pow(2, max_depth);

  double xi = a;
  for ( int i = 0; i < nbase; i++, xi += dx ) {
    xx.push_back(xi);
    yy.push_back(target(xi, data));
    zz.push_back(i > 0);
  }
  compute_spline();

  bool flag = true;
  while ( flag )
    flag = refine();
  return !flag; // will always be true, as failure causes error...
}

// Takes the contents of the 'xx' and 'yy' arrays and computes a
// (non-adaptive) spline with them.
void AdaptiveSpline::compute_spline() {
  const int n = xx.size();
  x.resize(n);
  y.resize(n);
  std::copy(xx.begin(), xx.end(), x.begin());
  std::copy(yy.begin(), yy.end(), y.begin());
  init_self();
}

// Refine the spline by adding points in all intervals (xx[i-1],
// xx[i]) where z[i] is true.
bool AdaptiveSpline::refine() {
  dx /= 2;

  if ( dx < dxmin )
    Rf_error("Spline as refined as currently possible");

  bool flag = false;
  
  std::list<double>::iterator xi = xx.begin(), yi = yy.begin();
  std::list<bool>::iterator zi = zz.begin();
  for ( ; xi != xx.end(); ++xi, ++yi, ++zi ) {
    if ( *zi ) {
      const double x_mid = *xi - dx;
      const double y_mid = target(x_mid, data);
      const double p_mid = eval(x_mid);

      // Always insert the new points.
      xx.insert(xi, x_mid);
      yy.insert(yi, y_mid);

      // Flag for refinement based on
      const bool flag_mid = !check_err(y_mid, p_mid);
      // If error was OK/not OK (flag_mid is true/false), say that
      // this interval is OK/not OK...
      *zi = flag_mid;
      // ...and that that the interval implied by the mid point is
      // OK/not OK.
      zz.insert(zi, flag_mid);

      flag = flag || flag_mid;
    }
  }

  // Recompute spline to use new points added during refinement.
  compute_spline();

  return flag;
}

// Determines if difference between predicted and true values falls
// within error bounds.
bool AdaptiveSpline::check_err(double y_true, double y_pred) const {
  const double err_abs = fabs(y_true - y_pred);
  const double err_rel = fabs(1 - y_pred / y_true);
  return err_abs < atol || err_rel < rtol;
}

// Get everything cleaned up.
void AdaptiveSpline::reset() {
  Spline::reset();
  zz.clear();
  xx.clear();
  yy.clear();
}

void AdaptiveSpline::check_bounds() {
  if ( a >= b )
    Rf_error("Impossible bounds");
  if ( !is_finite(a) || !is_finite(b) )
    Rf_error("Infinite bounds");
}
