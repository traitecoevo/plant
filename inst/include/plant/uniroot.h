// -*-c++-*-
#ifndef PLANT_PLANT_UNIROOT_H_
#define PLANT_PLANT_UNIROOT_H_

// Really simple wrapper around Boost's 1d root finding with bisection
// method.

#include <boost/math/tools/roots.hpp>
#include <plant/util.h>

namespace plant {
namespace util {

namespace internals {
struct uniroot_tol {
  uniroot_tol(double atol_, double rtol_) : atol(atol_), rtol(rtol_) {}
  bool operator()(double a, double b) {
    return std::abs(a - b) < atol + rtol * std::min(std::abs(a), std::abs(b));
  }
  double atol;
  double rtol;
};
}

// Wrapper around boost's root finder as a black-box function.
template <typename Function>
double uniroot(Function f, double min, double max, double tol,
               size_t max_iterations) {
  using boost::math::tools::bisect;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root = bisect(f, min, max,
                                          internals::uniroot_tol(tol, tol),
                                          it);
  if (it > static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop("Exceeded max_iterations");
  }
  return (root.first + root.second) / 2.0;
}

// Faster root finder for SMOOTH, monotonic functions (TOMS748 / Brent-like).
//
// Drop-in replacement for uniroot() with the same [min, max] bracketing
// contract (f(min) and f(max) must have opposite signs). For smooth functions
// TOMS748 converges super-linearly (~5-8 evals) versus bisection's ~one bit per
// iteration, so it is attractive for deeply nested, expensive-per-eval solvers.
//
// WARNING (empirical, root_water_uptake branch): substituting this for
// util::uniroot in the leaf hydraulic solvers (find_root_psi / psi_stem_to_ci)
// DESTABILISED the coupled soil-water ODE - it produced NaN soil potentials and
// ran *slower* overall. The continuum functions there are not smooth enough
// (vulnerability-curve clamps, splines with extrapolation disabled, near-flat
// regions), so the interpolation steps stall or probe bad points. Bisection's
// robustness is load-bearing in that code path. Use this only where the target
// function is known to be smooth and well-behaved across the whole bracket.
template <typename Function>
double uniroot_smooth(Function f, double min, double max, double tol,
                      size_t max_iterations) {
  using boost::math::tools::toms748_solve;
  boost::uintmax_t it = max_iterations;
  std::pair<double, double> root = toms748_solve(
      f, min, max, internals::uniroot_tol(tol, tol), it);
  if (it >= static_cast<boost::uintmax_t>(max_iterations)) {
    util::stop("Exceeded max_iterations");
  }
  return (root.first + root.second) / 2.0;
}

}
}

#endif
