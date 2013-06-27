#include "adaptive_spline.h"

namespace spline {

// Constructor that sets things to sensible defaults.  Setting the
// control parameters (via `set_control`) is optional, and the values
// here should be reasonable for many uses.
AdaptiveSpline::AdaptiveSpline(util::DFunctor *target)
  : target(target),
    atol(1e-6),
    rtol(1e-6),
    dx(NA_REAL),
    dxmin(NA_REAL),
    nbase(17),
    max_depth(16) {
}

// Set *all* the control parameters.
void AdaptiveSpline::set_control(double atol_, double rtol_, 
				 int nbase_, int max_depth_) {
  atol = atol_;
  rtol = rtol_;
  nbase = nbase_;
  max_depth = max_depth_;
}

// Evaluate the underlying function (double -> double).
double AdaptiveSpline::eval_target(double x) const {
  return (*target)(x);
}

// Adaptively refine a spline that spans from a to b so that the
// midpoints of evaluated points have sufficiently low error.
Spline AdaptiveSpline::construct_spline(double a, double b) {
  check_bounds(a, b);
  dx = (b - a) / (nbase - 1);
  dxmin = dx / pow(2, max_depth);

  xx.clear();
  yy.clear();
  zz.clear();

  double xi = a;
  for ( int i = 0; i < nbase; i++, xi += dx ) {
    xx.push_back(xi);
    yy.push_back(eval_target(xi));
    zz.push_back(i > 0);
  }

  compute_spline();

  bool flag = true;
  while ( flag )
    flag = refine();

  return spline;
}

// Given our current set of x/y points, construct the spline.  This is
// the only function that actually modifies the spline object.
void AdaptiveSpline::compute_spline() {
  std::vector<double> x(xx.begin(), xx.end());
  std::vector<double> y(yy.begin(), yy.end());
  spline.init(x, y);
}

// Refine the spline by adding points in all intervals (xx[i-1],
// xx[i]) where z[i] is true, and check to see if this new point has
// sufficiently good error that this interval needs no further
// refinement.
bool AdaptiveSpline::refine() {
  dx /= 2;

  if ( dx < dxmin )
    ::Rf_error("Spline as refined as currently possible");

  bool flag = false;
  
  std::list<double>::iterator xi = xx.begin(), yi = yy.begin();
  std::list<bool>::iterator zi = zz.begin();
  for ( ; xi != xx.end(); ++xi, ++yi, ++zi ) {
    if ( *zi ) {
      const double x_mid = *xi - dx;
      const double y_mid = eval_target(x_mid);
      const double p_mid = spline.eval(x_mid);

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

void AdaptiveSpline::check_bounds(double a, double b) {
  if ( a >= b )
    ::Rf_error("Impossible bounds");
  if ( !util::is_finite(a) || !util::is_finite(b) )
    ::Rf_error("Infinite bounds");
}

// Determines if difference between predicted and true values falls
// within error bounds.
bool AdaptiveSpline::check_err(double y_true, double y_pred) const {
  const double err_abs = fabs(y_true - y_pred);
  const double err_rel = fabs(1 - y_pred / y_true);
  return err_abs < atol || err_rel < rtol;
}


}
