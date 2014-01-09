#include "spline.h"

namespace spline {

// On contruction, make sure that `interp_accel` and `spline` are
// NULL, as we'll use these values to determine if memory needs
// clearing.
Spline::Spline()
  : interpolation_type(select_interpolation_type(false, false)),
    interp(NULL),
    interp_accel(NULL) {
}

Spline::Spline(bool is_akima, bool is_linear)
  : interpolation_type(select_interpolation_type(is_akima, is_linear)),
    interp(NULL),
    interp_accel(NULL) {
}

// On cleanup, remove memory allocated for the GSL support.
Spline::~Spline() {
  gsl_free_interp();
  gsl_free_interp_accel();
}

// Copy constructor -- needed because we need to reallocate the memory
// used by `interp_accel` and `spline`.  This is a bit annoying,
// because it involves recomputing the spline parameters.  A deep copy
// of `interp_accel` and `spline` would avoid this, but require more
// careful inspection of what GSL is storing (and dependence on the
// underlying implementation).
Spline::Spline(const Spline& obj)
  : x(obj.x),
    y(obj.y),
    interpolation_type(obj.interpolation_type),
    interp(NULL),
    interp_accel(NULL) {
  if (x.size() > 0)
    init_self();
}

Spline& Spline::operator=(Spline other) {
  using std::swap;
  swap(x,                  other.x);
  swap(y,                  other.y);
  swap(interpolation_type, other.interpolation_type);
  swap(interp,             other.interp);
  swap(interp_accel,       other.interp_accel);
  return *this;
}

std::string Spline::type() const {
  check_initialised();
  std::string ret(gsl_interp_name(interp));
  return ret;
}

// Build a spline out of the vectors 'x' and 'y'.
void Spline::init(const std::vector<double>& x_,
		  const std::vector<double>& y_) {
  util::check_length(y_.size(), x_.size());
  x = x_;
  y = y_;
  init_self();
}

// Compute the spline from the points contained in 'x' and 'y'.
void Spline::init_self() {
  const size_t n = x.size();

  if (interp == NULL) {
    interp = gsl_interp_alloc(interpolation_type, n);
  } else if (interp->size != n) {
    gsl_free_interp();
    interp = gsl_interp_alloc(interpolation_type, n);
  }
  gsl_interp_init(interp, ref(x), ref(y), n);

  gsl_free_interp_accel();
  interp_accel = gsl_interp_accel_alloc();
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
  gsl_free_interp();
  gsl_free_interp_accel();
}

// Compute the value of the spline at point `x=u`
double Spline::eval(double u) const {
  check_initialised();
  return gsl_interp_eval(interp, ref(x), ref(y), u, interp_accel);
}

double Spline::deriv(double u) const {
  check_initialised();
  return gsl_interp_eval_deriv(interp, ref(x), ref(y), u, interp_accel);
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
Rcpp::NumericMatrix Spline::r_get_xy() const {
  std::vector< std::vector<double> > xy;
  xy.push_back(x);
  xy.push_back(y);
  return util::to_rcpp_matrix(xy);
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
void Spline::gsl_free_interp() {
  if (interp != NULL) {
    gsl_interp_free(interp);
    interp = NULL;
  }
}

void Spline::gsl_free_interp_accel() {
  if (interp_accel != NULL) {
    gsl_interp_accel_free(interp_accel);
    interp_accel = NULL;
  }
}

void Spline::check_initialised() const {
  if (interp == NULL)
    Rcpp::stop("Interpolator not initialised -- cannot evaluate");
}

const gsl_interp_type*
Spline::select_interpolation_type(bool is_akima, bool is_linear) {
  if (is_linear) {
    return gsl_interp_linear;
  } else if (is_akima) {
    return gsl_interp_akima;
  } else {
    return gsl_interp_cspline;
  }
  return NULL;
}

}
