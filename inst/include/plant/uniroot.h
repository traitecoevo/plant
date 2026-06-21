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
// HISTORY (empirical, root_water_uptake branch, #486): an early *blanket* swap
// of this for util::uniroot across BOTH nested leaf hydraulic root-finds at once
// destabilised the coupled soil-water ODE (NaN soil potentials, slower overall),
// which was first read as "the hydraulic path is too non-smooth for a
// superlinear solver". That conclusion was too broad. Re-examined target by
// target, both leaf solvers are in fact smooth and strictly monotone over the
// brackets they are actually handed, and both now use this solver at their
// existing tolerances (same root, fewer evals):
//   * psi_stem_to_ci (Phase 6): A_colim demand minus the linear gc supply over
//     (gamma*, ca]; ~29 -> ~9 evals at 1e-7.
//   * find_root_psi (Phase 8): the soil->collar continuity residual over
//     [-psi_crit, wettest_soil_layer]; ~15-16 -> ~6-8 evals at 1e-4. Its
//     brackets are guaranteed opposite-sign/finite by find_root_collar_psi's
//     early-exits.
// The genuine non-smoothness (vulnerability-curve clamps, the root vulnerability
// spline extrapolating negative beyond its domain, near-flat regions) lives in
// E_from_Soil_to_Root_Collar itself, NOT in the root-finders, and would break
// bisection too. So: use this where the target is smooth and well-behaved across
// the whole bracket -- which the leaf solvers are, on their operating brackets.
// One gotcha vs bisect: this validates its bracket and THROWS on non-finite or
// same-sign endpoints where boost::bisect returned NaN silently (guard upstream
// if a finite-but-degenerate bracket can occur; see psi_stem_to_ci's NA guard).
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
