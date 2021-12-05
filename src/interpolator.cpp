#include <plant/interpolator.h>
#include <plant/util.h>
#include <plant/util_post_rcpp.h> // to_rcpp_matrix
#include <Rcpp.h>

namespace plant {
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
  // https://stackoverflow.com/questions/17769114/stdis-sorted-and-strictly-less-comparison
  if (not std::is_sorted(x.begin(), x.end(), std::less_equal<double>())) {
    util::stop("spline control points must be unique and in ascending order");
  }
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

// adds point in sorted position (slower than above)
void Interpolator::add_point_sorted(double xi, double yi) {
  auto x_upper = std::upper_bound(x.begin(), x.end(), xi); // find smallest number larger than xi
  x.insert(x_upper, xi); // add xi below that number
  auto y_upper = std::upper_bound(y.begin(), y.end(), yi);
  y.insert(y_upper, yi);
}

// Remove all the contents, being ready to be refilled.
void Interpolator::clear() {
  x.clear();
  y.clear();
  active = false;
}

// Compute the value of the interpolated function at point `x=u`
double Interpolator::eval(double u) const {
  check_active(); // slower to check this every time
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
  auto ret = std::vector<double>();
  ret.reserve(u.size()); // fast to do this once rather than multiple times with push_back
  for (auto const& x : u) {
    ret.push_back(eval(x));
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
