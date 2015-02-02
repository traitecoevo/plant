#include <tree2/interpolator.h>
#include <tree2/util.h>
#include <tree2/util_post_rcpp.h> // to_rcpp_matrix
#include <Rcpp.h>

namespace tree2 {
namespace interpolator {

// Build a interpolator out of the vectors 'x' and 'y'.
void Interpolator_TK::init(const std::vector<double>& x_,
			const std::vector<double>& y_) {
  util::check_length(y_.size(), x_.size());
  x = x_;
  y = y_;
  initialise();
}

// Compute the interpolated function from the points contained in 'x'
// and 'y'.
void Interpolator_TK::initialise() {
  tk_spline.set_points(x, y);
}

// Support for adding points in turn (assumes monotonic increasing in
// 'x', unchecked).
void Interpolator_TK::add_point(double xi, double yi) {
  x.push_back(xi);
  y.push_back(yi);
}

// Remove all the contents, being ready to be refilled.
void Interpolator_TK::clear() {
  x.clear();
  y.clear();
}

// Compute the value of the interpolated function at point `x=u`
double Interpolator_TK::eval(double u) const {
  return tk_spline(u);
}

double Interpolator_TK::deriv(double /* u */) const {
  util::stop("Derivatives not implemented");
  return NA_REAL;
}

// Return the number of (x,y) pairs contained in the Interpolator_TK.
size_t Interpolator_TK::size() const {
  return x.size();
}

// These are chosen so that if a Interpolator_TK is empty, functions
// looking to see if they will fall outside of the covered range will
// always find they do.  This is the same principle as R's
// range(numeric(0)) -> c(Inf, -Inf)
double Interpolator_TK::min() const {
  return size() > 0 ? x.front() : R_PosInf;
}
double Interpolator_TK::max() const {
  return size() > 0 ? x.back() : R_NegInf;
}

std::vector<double> Interpolator_TK::get_x() const {
  return x;
}

std::vector<double> Interpolator_TK::get_y() const {
  return y;
}

// Get the (x,y) pairs in the Interpolator_TK as a two-column matrix
SEXP Interpolator_TK::r_get_xy() const {
  std::vector< std::vector<double> > xy;
  xy.push_back(x);
  xy.push_back(y);
  return Rcpp::wrap(util::to_rcpp_matrix(xy));
}

// Compute the value of the interpolated function at a vector of
// points `x=u`, returning a vector of the same length.
std::vector<double> Interpolator_TK::r_eval(std::vector<double> u) const {
  const size_t n = u.size();
  std::vector<double> ret(n);
  for (size_t i = 0; i < n; ++i) {
    ret[i] = eval(u[i]);
  }
  return ret;
}

std::vector<double> Interpolator_TK::r_deriv(std::vector<double> u) const {
  const size_t n = u.size();
  std::vector<double> ret(n);
  for (size_t i = 0; i < n; ++i) {
    ret[i] = deriv(u[i]);
  }
  return ret;
}

}
}
