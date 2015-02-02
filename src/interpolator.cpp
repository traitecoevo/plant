#include <tree2/interpolator.h>
#include <tree2/util.h>
#include <tree2/util_post_rcpp.h> // to_rcpp_matrix
#include <Rcpp.h>

namespace tree2 {
namespace interpolator {

// Build a interpolator out of the vectors 'x' and 'y'.
void Interpolator::init(const std::vector<double>& x_,
			const std::vector<double>& y_) {
  util::check_length(y_.size(), x_.size());
  if (x_.size() < 3) {
    util::stop("insufficient number of points");
  }
  x = x_;
  y = y_;
  initialise();
}

// Compute the interpolated function from the points contained in 'x'
// and 'y'.
void Interpolator::initialise() {
  if (x.size() > 0) {
    tk_spline.set_points(x, y);
    active = true;
  }
}

// Support for adding points in turn (assumes monotonic increasing in
// 'x', unchecked).
void Interpolator::add_point(double xi, double yi) {
  x.push_back(xi);
  y.push_back(yi);
}

// Remove all the contents, being ready to be refilled.
void Interpolator::clear() {
  x.clear();
  y.clear();
  active = false;
}

// Compute the value of the interpolated function at point `x=u`
double Interpolator::eval(double u) const {
  return tk_spline(u);
}

// Return the number of (x,y) pairs contained in the Interpolator.
size_t Interpolator::size() const {
  return x.size();
}

// These are chosen so that if a Interpolator is empty, functions
// looking to see if they will fall outside of the covered range will
// always find they do.  This is the same principle as R's
// range(numeric(0)) -> c(Inf, -Inf)
double Interpolator::min() const {
  return size() > 0 ? x.front() : R_PosInf;
}
double Interpolator::max() const {
  return size() > 0 ? x.back() : R_NegInf;
}

std::vector<double> Interpolator::get_x() const {
  return x;
}

std::vector<double> Interpolator::get_y() const {
  return y;
}

// Get the (x,y) pairs in the Interpolator as a two-column matrix
SEXP Interpolator::r_get_xy() const {
  std::vector< std::vector<double> > xy;
  xy.push_back(x);
  xy.push_back(y);
  return Rcpp::wrap(util::to_rcpp_matrix(xy));
}

// Compute the value of the interpolated function at a vector of
// points `x=u`, returning a vector of the same length.
std::vector<double> Interpolator::r_eval(std::vector<double> u) const {
  check_active();
  const size_t n = u.size();
  std::vector<double> ret(n);
  for (size_t i = 0; i < n; ++i) {
    ret[i] = eval(u[i]);
  }
  return ret;
}

void Interpolator::check_active() const {
  if (!active) {
    util::stop("Interpolator not initialised -- cannot evaluate");
  }
}

}
}
