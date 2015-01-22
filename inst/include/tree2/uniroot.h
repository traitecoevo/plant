// -*-c++-*-
#ifndef TREE_UNIROOT_H_
#define TREE_UNIROOT_H_

// Really simple wrapper around Boost's 1d root finding with bisection
// method.

#include <boost/math/tools/roots.hpp>

namespace tree2 {
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
  std::pair<double, double> root = bisect(f, min, max,
                                          internals::uniroot_tol(tol, tol),
                                          max_iterations);
  // TODO: probably should check here that we didn't exhaust the
  // number of steps we may take (i.e. max_iterations was not too bad).
  return (root.first + root.second) / 2.0;
}

}
}

#endif
