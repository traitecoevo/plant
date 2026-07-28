#include <plant/adaptive_interpolator.h>
#include <Rcpp.h>
#include <plant/util_post_rcpp.h>
#include <cmath>

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

// Find the interval refinement is stuck on. zz[i] flags the interval
// (xx[i-1], xx[i]) as still failing the error test; every one of them is about
// 2*dx wide at this point, so the informative one is the one across which the
// target moves furthest -- that is the feature refinement cannot resolve.
StallReport AdaptiveInterpolator::stall_report() const {
  StallReport ret = {xx.front(), xx.back(), yy.front(), yy.back(), 0};
  double jump_max = -1.0;

  std::list<double>::const_iterator xi = xx.begin(), yi = yy.begin();
  std::list<bool>::const_iterator zi = zz.begin();
  double x_prev = *xi, y_prev = *yi;

  for (++xi, ++yi, ++zi; xi != xx.end(); ++xi, ++yi, ++zi) {
    if (*zi) {
      ret.n_unresolved++;
      const double jump = fabs(*yi - y_prev);
      if (jump > jump_max) {
        jump_max = jump;
        ret.x_lo = x_prev;
        ret.x_hi = *xi;
        ret.y_lo = y_prev;
        ret.y_hi = *yi;
      }
    }
    x_prev = *xi;
    y_prev = *yi;
  }

  return ret;
}

// Report a resolution limit, naming the interval it was hit on.
//
// The bare "as refined as currently possible" message says only that some
// feature is too narrow, which leaves the reader to find it by hand -- an
// afternoon's work in the case that prompted this (#571: a dryland patch whose
// cohorts had all stalled at the introduction height, putting a step in the
// light profile). The location turns that into a one-line diagnosis, and callers
// with model context add to it by catching refinement_failure.
[[noreturn]] void AdaptiveInterpolator::stop_unresolvable() const {
  const StallReport r = stall_report();

  std::string msg =
      "Interpolated function as refined as currently possible: spacing " +
      util::format_double(dx) + " is below the limit " +
      util::format_double(dxmin) + " set by max_depth=" +
      util::to_string(static_cast<int>(max_depth)) +
      ", and the target still misses atol=" + util::format_double(atol) +
      " / rtol=" + util::format_double(rtol) +
      ". The target has a feature narrower than that spacing";

  if (r.n_unresolved > 0) {
    msg += ", at x = " + util::format_double(r.x_lo) + ": across [" +
           util::format_double(r.x_lo) + ", " + util::format_double(r.x_hi) +
           "] the target jumps from " + util::format_double(r.y_lo) + " to " +
           util::format_double(r.y_hi) + " (" +
           util::to_string(r.n_unresolved) +
           " interval(s) still unresolved, over the domain [" +
           util::format_double(xx.front()) + ", " +
           util::format_double(xx.back()) + "])";
  }
  msg += ".";

  throw refinement_failure(msg, r);
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
