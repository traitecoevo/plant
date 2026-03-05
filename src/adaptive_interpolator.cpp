#include <plant/adaptive_interpolator.h>
#include <Rcpp.h>
#include <plant/util_post_rcpp.h>

namespace plant {
namespace interpolator {

// Given our current set of x/y points, construct the interpolated
// fucntion.  This is the only function that actually modifies the
// interpolator object.
void AdaptiveInterpolator::update_spline() {
  std::vector<double> x(xx.begin(), xx.end());
  std::vector<double> y(yy.begin(), yy.end());
  interpolator.init(x, y);
}

void AdaptiveInterpolator::check_bounds(double a, double b) {
  if (a >= b) {
    Rcpp::stop("Impossible bounds");
  }
  if (!util::is_finite(a) || !util::is_finite(b)) {
    Rcpp::stop("Infinite bounds");
  }
}

// Determines if difference between predicted and true values falls
// within error bounds.
bool AdaptiveInterpolator::check_err(double y_true, double y_pred) const {
  const double err_abs = fabs(y_true - y_pred);
  const double err_rel = fabs(1 - y_pred / y_true);
  return err_abs < atol || err_rel < rtol;
}

}
}
