// -*-c++-*-
#ifndef TREE_GRADIENT_H_
#define TREE_GRADIENT_H_

#include <vector>
#include <cmath>
#include <cstddef>
#include <limits>

namespace tree2 {
namespace util {

// A. One-shot

// 1. Forward difference:
template <typename Function>
double gradient_fd_forward(Function f, double x, double dx) {
  return gradient_fd_forward(f, x, dx, f(x));
}
template <typename Function>
double gradient_fd_forward(Function f, double x, double dx, double fx) {
  return (f(x + dx) - fx) / dx;
}

// 2. Backward difference (just wraps around forward difference with
// the direction flipped)
template <typename Function>
double gradient_fd_backward(Function f, double x, double dx) {
  return gradient_fd_forward(f, x, -dx);
}

template <typename Function>
double gradient_fd_backward(Function f, double x, double dx, double fx) {
  return gradient_fd_forward(f, x, -dx, fx);
}

// 3. Centre (can't use f(x))
template <typename Function>
double gradient_fd_centre(Function f, double x, double dx) {
  const double dx2 = dx / 2;
  return (f(x + dx2) - f(x - dx2)) / dx;
}

// 4. Wrapper:
template <typename Function>
double gradient_fd(Function f, double x, double dx, int direction) {
  if (direction < 0) {
    return gradient_fd_backward(f, x, dx);
  } else if (direction == 0) {
    return gradient_fd_centre(f, x, dx);
  } else {
    return gradient_fd_forward(f, x, dx);
  }
}

template <typename Function>
double gradient_fd(Function f, double x, double dx, double fx, int direction) {
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
template <typename Function>
double gradient_richardson(Function f, double x, double d, size_t r) {
  const size_t v = 2; // this value is required by scheme (above)
  const double zero_tol = sqrt(std::numeric_limits<double>::epsilon())/7e-7;

  // Initial offset (see above).
  double h = std::abs(d * x) + d * (std::abs(x) < zero_tol);

  std::vector<double> a;
  for (size_t i = 0; i < r; i++, h /= v) {
    a.push_back((f(x + h) - f(x - h))/(2*h));
  }

  for (size_t m = 1; m < r; ++m) {
    const double four_m = pow(4.0, m);
    std::vector<double> a_next;
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
