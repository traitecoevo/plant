// -*-c++-*-
#ifndef PLANT_PLANT_GRADIENT_H_
#define PLANT_PLANT_GRADIENT_H_

#include <vector>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace plant {
namespace util {

// The scalar is deduced from the point the derivative is taken at, so a double
// caller instantiates at double and gets the arithmetic that was here. The step
// stays double: it is a spacing, not a position, and nothing differentiates it.

// The integrand has to return the scalar of the point it is evaluated at. An
// integrand returning double at an active point takes the value of an active
// quantity, so every quotient below is exactly zero while the result's type
// stays active and nothing is raised.
template <typename Function, typename S>
concept integrand_of =
  std::same_as<std::invoke_result_t<Function&, const S&>, S>;

// A. One-shot

// 1. Forward difference:
template <typename Function, typename S>
  requires integrand_of<Function, S>
S gradient_fd_forward(Function f, const S& x, double dx, const S& fx) {
  return (f(x + dx) - fx) / dx;
}
template <typename Function, typename S>
  requires integrand_of<Function, S>
S gradient_fd_forward(Function f, const S& x, double dx) {
  return gradient_fd_forward(f, x, dx, S(f(x)));
}

// 2. Backward difference (just wraps around forward difference with
// the direction flipped)
template <typename Function, typename S>
  requires integrand_of<Function, S>
S gradient_fd_backward(Function f, const S& x, double dx) {
  return gradient_fd_forward(f, x, -dx);
}

template <typename Function, typename S>
  requires integrand_of<Function, S>
S gradient_fd_backward(Function f, const S& x, double dx, const S& fx) {
  return gradient_fd_forward(f, x, -dx, fx);
}

// 3. Centre (can't use f(x))
template <typename Function, typename S>
  requires integrand_of<Function, S>
S gradient_fd_centre(Function f, const S& x, double dx) {
  const double dx2 = dx / 2;
  return (f(x + dx2) - f(x - dx2)) / dx;
}

// 4. Wrapper:
template <typename Function, typename S>
  requires integrand_of<Function, S>
S gradient_fd(Function f, const S& x, double dx, int direction) {
  if (direction < 0) {
    return gradient_fd_backward(f, x, dx);
  } else if (direction == 0) {
    return gradient_fd_centre(f, x, dx);
  } else {
    return gradient_fd_forward(f, x, dx);
  }
}

template <typename Function, typename S>
  requires integrand_of<Function, S>
S gradient_fd(Function f, const S& x, double dx, const S& fx, int direction) {
  if (direction < 0) {
    return gradient_fd_backward(f, x, dx, fx);
  } else if (direction == 0) {
    return gradient_fd_centre(f, x, dx); // no fx
  } else {
    return gradient_fd_forward(f, x, dx, fx);
  }
}

// B. Multi-shot

// Richardson extrapolation of centred difference -- the most accurate
// but the slowest to compute.
//
// Based on code in R's numDeriv::grad (actually in grad.default).
//
// First order derivatives are stored in the vector a[r], for r rounds
// of improvement.
//
// We start with deviation from x of d * x, unless x is almost zero
// (determined by being smaller than zero_tol) in which case we use
// d as the absolute deviation.

//------------------------------------------------------------------------
//   Applying Richardson Extrapolation to improve the accuracy of
//   the first and second order derivatives. The algorithm as follows:
//
//   --  For each column of the derivative matrix a,
//        say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
//        new sequence of approximations B1, B2, ..., Br used the formula
//
//           B(i) =( A(i+1)*4^m - A(i) ) / (4^m - 1) ,  i=1,2,...,r-m
//
//              N.B. This formula assumes v=2.
//
//   -- Initially m is taken as 1  and then the process is repeated
//       restarting with the latest improved values and increasing the
//       value of m by one each until m equals r-1
//
//-------------------------------------------------------------------------
template <typename Function, typename S>
  requires integrand_of<Function, S>
S gradient_richardson(Function f, const S& x, double d, size_t r) {
  const size_t v = 2; // this value is required by scheme (above)
  const double zero_tol = sqrt(std::numeric_limits<double>::epsilon())/7e-7;

  // Unqualified so an active scalar reaches its own abs by name, and std::abs
  // still answers for double.
  using std::abs;
  // Initial offset (see above).
  S h = abs(d * x) + d * (abs(x) < zero_tol);

  std::vector<S> a;
  for (size_t i = 0; i < r; i++, h /= v) {
    a.push_back((f(x + h) - f(x - h))/(2*h));
  }

  for (size_t m = 1; m < r; ++m) {
    const double four_m = pow(4.0, m);
    std::vector<S> a_next;
    for (size_t i = 0; i < r - m; ++i) {
      a_next.push_back((a[i+1]*four_m - a[i])/(four_m - 1));
    }
    a = a_next;
  }

  return a.front();
}

}
}

#endif
