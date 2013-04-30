#include "spline.h"

namespace spline {

// On contruction, make sure that `acc` and `spline` are NULL, as
// we'll use these values to determine if memory needs clearing.
Spline::Spline() {
  acc    = NULL;
  spline = NULL;
}

// On cleanup, remove memory allocated for the GSL support.
Spline::~Spline() {
  gsl_free_spline();
  gsl_free_acc();
}

// Copy constructor -- needed because we need to reallocate the memory
// used by `acc` and `spline`.  This is a bit annoying, because it
// involves recomputing the spline parameters.  A deep copy of `acc`
// and `spline` would avoid this, but require more careful inspection
// of what GSL is storing (and dependence on the underlying
// implementation).
Spline::Spline(const Spline& obj) : x(obj.x), y(obj.y) {
  acc = NULL;
  spline = NULL;
  if ( x.size() > 0 )
    init_self();
}

// Build a spline out of the vectors 'x' and 'y'.  
// TODO: This might be nicer if we did this with iterators, or with a
// template?
void Spline::init(std::vector<double> x_, std::vector<double> y_) {
  x      = x_;
  y      = y_;
  init_self();
}

// Compute the spline from the points contained in 'x' and 'y'.
void Spline::init_self() {
  const size_t n = x.size();
  gsl_free_spline();
  gsl_free_acc();
  acc    = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, &x.front(), &y.front(), n);  
}

// Support for adding points in turn (assumes monotonic increasing in
// 'x', unchecked).
void Spline::add_point(double xi, double yi) {
  x.push_back(xi);
  y.push_back(yi);
}

// Remove all the contents, being ready to be refilled.
void Spline::reset() {
  x.clear();
  y.clear();
  gsl_free_spline();
}

// Compute the value of the spline at point `x=u`
// TODO: this will crash if called before `init`.
double Spline::eval(double u) const {
  return gsl_spline_eval(spline, u, acc);
}

// Return the number of (x,y) pairs contained in the spline.
size_t Spline::size() const {
  return x.size();
}

Rcpp::NumericVector Spline::r_get_x() const {
  Rcpp::NumericVector ret(x.begin(), x.end());
  return ret;
}

Rcpp::NumericVector Spline::r_get_y() const {
  Rcpp::NumericVector ret(y.begin(), y.end());
  return ret;
}

// Get the (x,y) pairs in the spline as a two-column matrix
// NOTE: line 2 `ret(n, 2)` causes an error, but I don't see how to
// cast my way out of it.
Rcpp::NumericMatrix Spline::r_get_xy() const {
  const size_t n = x.size();
  Rcpp::NumericMatrix ret(n, 2);

  for ( size_t i = 0; i < n; i++ ) {
    ret(i,0) = x[i];
    ret(i,1) = y[i];
  }
  
  return ret;
}

// Compute the value of the spline at a vector of points `x=u`,
// returning a vector of the same length.
std::vector<double> Spline::r_eval(std::vector<double> u) const {
  const size_t n = u.size();
  std::vector<double> ret(n);
  for ( size_t i = 0; i < n; i++ )
    ret[i] = eval(u[i]);
  return ret;
}

// Helper functions to cleanup GSL memory if it looks like it can be
// cleaned up.
void Spline::gsl_free_spline() {
  if ( spline != NULL ) {
    gsl_spline_free(spline);
    spline = NULL;
  }
}

void Spline::gsl_free_acc() {
  if ( acc != NULL ) {
    gsl_interp_accel_free(acc);
    acc = NULL;
  }
}

}
