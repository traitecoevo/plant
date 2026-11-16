// -*-c++-*-
#ifndef PLANT_PLANT_UNIROOT_H_
#define PLANT_PLANT_UNIROOT_H_

// A bracketed one-dimensional root-find: Boost's TOMS748 (Alefeld, Potra and
// Shi 1995), as phylloptim uses for its smooth residuals.

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

// More than a smooth target takes by orders of magnitude, so reaching it means
// the target is not one this should be used on.
constexpr size_t uniroot_max_iterations = 1000;

// The root of f in [min, max], which f must bracket (opposite signs at the
// ends), as the midpoint of a final bracket narrower than tol * (1 + |root|).
// TOMS748 keeps a bracket like bisection, so it cannot leave [min, max], and on
// a smooth target converges superlinearly where bisection gains one bit per
// step. Refuses where the bracket never reached the tolerance: Boost reports
// the iterations it used, which never exceed the cap, so counting them cannot
// tell a converged search from an exhausted one.
template <typename Function>
double uniroot(Function f, double min, double max, double tol) {
  using boost::math::tools::toms748_solve;
  internals::uniroot_tol converged(tol, tol);
  boost::uintmax_t it = uniroot_max_iterations;
  const std::pair<double, double> root =
      toms748_solve(f, min, max, converged, it);
  if (!(root.first == root.second || converged(root.first, root.second))) {
    util::stop("uniroot: the bracket [" + util::to_string(root.first) + ", " +
               util::to_string(root.second) + "] is still wider than the "
               "tolerance after " + util::to_string(it) + " iterations");
  }
  return (root.first + root.second) / 2.0;
}

}
}

#endif
