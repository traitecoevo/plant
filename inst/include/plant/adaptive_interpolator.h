// -*-c++-*-
#ifndef PLANT_PLANT_ADAPTIVE_INTERPOLATOR_H_
#define PLANT_PLANT_ADAPTIVE_INTERPOLATOR_H_

// This class creates an interpolator object, using adaptive refinement

#include <list>
#include <stdexcept>
#include <string>
#include <plant/util.h> // util::stop, util::seq_len
#include <odelia/interpolator.hpp>

namespace plant {
namespace interpolator {

// Where refinement gave up: the worst interval still failing the error test at
// max_depth, and how many others are in the same state. Reported because the
// location is the diagnosis -- refinement stalls on a feature of the target, and
// naming the x it stalled at points straight at whatever put the feature there.
struct StallReport {
  double x_lo, x_hi;   // the unresolved interval
  double y_lo, y_hi;   // target values at its ends
  size_t n_unresolved; // intervals still failing, including this one
};

// Thrown when adaptive refinement exhausts max_depth. Carries the stall location
// so a caller with model context can say what lives at that x before the error
// reaches the user: the interpolator only sees a function with a narrow feature,
// whereas Patch knows it is a light profile over a set of cohort heights.
class refinement_failure : public std::runtime_error {
public:
  refinement_failure(const std::string& msg, const StallReport& report_)
    : std::runtime_error(msg), report(report_) {}
  StallReport report;
};

class AdaptiveInterpolator {
public:
  AdaptiveInterpolator(double atol_, double rtol_,
                       size_t nbase_, size_t max_depth_)
  : atol(atol_),
    rtol(rtol_),
    nbase(nbase_),
    max_depth(max_depth_),
    dx(NA_REAL),
    dxmin(NA_REAL),
    interpolator() {
  }
  AdaptiveInterpolator()
  : atol(0.0),
    rtol(0.0),
    nbase(0.0),
    max_depth(0.0),
    dx(NA_REAL),
    dxmin(NA_REAL),
    interpolator() {
  }

  template <typename Function>
  odelia::interpolator::Interpolator construct(Function target, double a, double b);
private:
  void update_spline();
  template <typename Function>
  bool refine(Function target);

  bool check_err(double y_true, double y_pred) const;
  static void check_bounds(double a, double b);
  static void check_target_value(double x, double y);

  // Locate the interval refinement is stuck on, and throw refinement_failure
  // naming it. Non-template members, so they live in the .cpp.
  StallReport stall_report() const;
  [[noreturn]] void stop_unresolvable() const;

  // Control parameters:
  double atol, rtol;
  size_t nbase, max_depth;
  double dx, dxmin;

  // In contrast to the Interpolator's x and y, which are vectors, these are
  // lists so that we can easily add points in the middle of them.
  std::list<double> xx, yy;
  std::list<bool> zz;

  // Temporary interpolator object that we build.
  odelia::interpolator::Interpolator interpolator;
};

// Adaptively refine a interpolator that spans from a to b so that the
// midpoints of evaluated points have sufficiently low error.
template <typename Function>
odelia::interpolator::Interpolator AdaptiveInterpolator::construct(Function target,
                                                                   double a, double b)
{
  check_bounds(a, b);
  dx = (b - a) / (static_cast<int>(nbase) - 1);
  dxmin = dx / pow(2, max_depth);

  xx.clear();
  yy.clear();
  zz.clear();

  std::vector<double> tmp = util::seq_len(a, b, nbase);
  for (size_t i = 0; i < nbase; i++) {
    const double xi = tmp[i];
    const double yi = target(xi);
    check_target_value(xi, yi);
    xx.push_back(xi);
    yy.push_back(yi);
    zz.push_back(i > 0);
  }

  update_spline();

  bool flag = true;
  while (flag) {
    flag = refine(target);
  }

  return interpolator;
}

// Refine the interpolated function by adding points in all intervals
// (xx[i-1], xx[i]) where z[i] is true, and check to see if this new
// point has sufficiently good error that this interval needs no
// further refinement.
template <typename Function>
bool AdaptiveInterpolator::refine(Function target) {
  dx /= 2;

  if (dx < dxmin) {
    // Say what was actually exhausted, and where. The target is finite here
    // (construct() and the loop below both check), so this is a genuine
    // resolution limit: the function has a feature narrower than dxmin, i.e. is
    // effectively discontinuous at this tolerance.
    stop_unresolvable();
  }

  bool flag = false;

  std::list<double>::iterator xi = xx.begin(), yi = yy.begin();
  std::list<bool>::iterator zi = zz.begin();
  for (; xi != xx.end(); ++xi, ++yi, ++zi) {
    if (*zi) {
      const double x_mid = *xi - dx;
      const double y_mid = target(x_mid);
      check_target_value(x_mid, y_mid);
      const double p_mid = interpolator.eval(x_mid);

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

  // Recompute interpolator to use new points added during refinement.
  update_spline();

  return flag;
}

}
}

#endif
