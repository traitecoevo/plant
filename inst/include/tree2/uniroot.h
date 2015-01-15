// -*-c++-*-
#ifndef TREE_UNIROOT_H_
#define TREE_UNIROOT_H_

// Really simple wrapper around Boost's 1d root finding with bisection
// method.

#include <boost/math/tools/roots.hpp>

namespace util {

namespace internals {
// TODO: Need to get both relative and absolute tolerance here I
// think, but I forget how I usually do that.  Not a biggie really.
struct uniroot_tol {
  uniroot_tol(double tol_) : tol(tol_) {}
  double tol;
  bool operator()(double min, double max) {
    return std::abs(min - max) < tol;
  }
};
}

// Wrapper around boost's root finder as a black-box function.
template <typename Function>
double uniroot(Function f, double min, double max, double tol,
               size_t max_iterations) {
  using boost::math::tools::bisect;
  std::pair<double, double> root = bisect(f, min, max,
                                          internals::uniroot_tol(tol),
                                          max_iterations);
  // TODO: probably should check here that we didn't exhaust the
  // number of steps we may take (i.e. max_iterations was not too bad).
  return (root.first + root.second) / 2.0;
}

}

#endif
