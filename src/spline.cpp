#include "spline.h"

namespace spline {

// On contruction, make sure that `acc` and `spline` are NULL, as
// we'll use these values to determine if memory needs clearing.
Spline::Spline()
  : acc(NULL),
    spline(NULL),
    is_akima(false) {
}

Spline::Spline(bool is_akima_)
  : acc(NULL),
    spline(NULL),
    is_akima(is_akima_) {
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
Spline::Spline(const Spline& obj)
  : x(obj.x),
    y(obj.y),
    acc(NULL),
    spline(NULL),
    is_akima(obj.is_akima) {
  if (x.size() > 0)
    init_self();
}

Spline& Spline::operator=(Spline other) {
  using std::swap;
  swap(x,        other.x);
  swap(y,        other.y);
  swap(acc,      other.acc);
  swap(spline,   other.spline);
  swap(is_akima, other.is_akima);
  return *this;
}

// Build a spline out of the vectors 'x' and 'y'.
void Spline::init(std::vector<double> x_, std::vector<double> y_) {
  x = x_;
  y = y_;
  init_self();
}

// Compute the spline from the points contained in 'x' and 'y'.
void Spline::init_self() {
  const size_t n = x.size();
  gsl_free_spline();
  gsl_free_acc();
  acc = gsl_interp_accel_alloc();
  if (is_akima) {
    spline = gsl_spline_alloc(gsl_interp_akima, n);
  } else {
    spline = gsl_spline_alloc(gsl_interp_cspline, n);
  }
  gsl_spline_init(spline, &x.front(), &y.front(), n);  
}

// Support for adding points in turn (assumes monotonic increasing in
// 'x', unchecked).
void Spline::add_point(double xi, double yi) {
  x.push_back(xi);
  y.push_back(yi);
}

// Remove all the contents, being ready to be refilled.
void Spline::clear() {
  x.clear();
  y.clear();
  gsl_free_spline();
}

// Compute the value of the spline at point `x=u`
double Spline::eval(double u) const {
  if (spline == NULL)
    ::Rf_error("Spline not initialised -- cannot evaluate");
  return gsl_spline_eval(spline, u, acc);
}

double Spline::deriv(double u) const {
  if (spline == NULL)
    ::Rf_error("Spline not initialised -- cannot evaluate");
  return gsl_spline_eval_deriv(spline, u, acc);
}

// Return the number of (x,y) pairs contained in the spline.
size_t Spline::size() const {
  return x.size();
}

// These are chosen so that if a spline is empty, functions looking to
// see if they will fall outside of the covered range will always find
// they do.  This is the same principle as R's 
//   range(numeric(0)) -> c(Inf, -Inf)
double Spline::min() const {
  return size() > 0 ? x.front() : R_PosInf;
}
double Spline::max() const {
  return size() > 0 ? x.back() : R_NegInf;
}

std::vector<double> Spline::get_x() const {
  return x;
}

std::vector<double> Spline::get_y() const {
  return y;
}

// Get the (x,y) pairs in the spline as a two-column matrix
// NOTE: line 2 `ret(n, 2)` causes an error, but I don't see how to
// cast my way out of it.
Rcpp::NumericMatrix Spline::r_get_xy() const {
  const size_t n = x.size();
  Rcpp::NumericMatrix ret(static_cast<int>(n), 2);

  for (size_t i = 0; i < n; ++i) {
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
  for (size_t i = 0; i < n; ++i)
    ret[i] = eval(u[i]);
  return ret;
}

std::vector<double> Spline::r_deriv(std::vector<double> u) const {
  const size_t n = u.size();
  std::vector<double> ret(n);
  for (size_t i = 0; i < n; ++i)
    ret[i] = deriv(u[i]);
  return ret;
}

// Helper functions to cleanup GSL memory if it looks like it can be
// cleaned up.
void Spline::gsl_free_spline() {
  if (spline != NULL) {
    gsl_spline_free(spline);
    spline = NULL;
  }
}

void Spline::gsl_free_acc() {
  if (acc != NULL) {
    gsl_interp_accel_free(acc);
    acc = NULL;
  }
}

}
