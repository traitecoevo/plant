#include <plant/adaptive_interpolator.h>
#include <Rcpp.h>
#include <plant/util_post_rcpp.h>

// Test hook: adaptively interpolate an R function, and return the knots chosen.
// The production caller (ResourceSpline, for the light profile) passes a C++
// lambda closed over a Patch, so there is no way to drive the refiner from R
// with a target of one's own choosing without this. Same pattern as
// test_gradient_fd1 in src/gradient.cpp.
//
// [[Rcpp::export]]
Rcpp::NumericVector test_adaptive_interpolator(Rcpp::Function f, double a,
                                               double b, double atol = 1e-6,
                                               double rtol = 1e-6,
                                               int nbase = 17,
                                               int max_depth = 16) {
  plant::interpolator::AdaptiveInterpolator ai(
      atol, rtol, static_cast<size_t>(nbase), static_cast<size_t>(max_depth));
  auto out = ai.construct(plant::util::RFunctionWrapper(f), a, b);
  return Rcpp::wrap(out.get_x());
}

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

// Refuse to refine against a non-finite target value.
//
// check_err() below compares NaN, and every comparison against NaN is false, so
// a single NaN or Inf from the target makes its interval permanently
// unacceptable: refinement halves the spacing until max_depth is exhausted and
// then reports "as refined as currently possible" — blaming resolution for what
// is really a non-finite value somewhere in the function being interpolated.
// That message has cost real debugging time (a TF24 patch whose competition
// profile went non-finite at low soil moisture read as an interpolation
// tolerance problem). Fail here instead, naming the point and the value.
void AdaptiveInterpolator::check_target_value(double x, double y) {
  if (!util::is_finite(y)) {
    Rcpp::stop("Adaptive interpolation target is non-finite (" +
               util::to_string(y) + ") at x = " + util::to_string(x) +
               "; refinement cannot fix a non-finite value.");
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
