// -*-c++-*-
#ifndef PLANT_PLANT_UTIL_H_
#define PLANT_PLANT_UTIL_H_

#include <stddef.h> // size_t
#include <RcppCommon.h> // as/wrap/SEXP
#include <R_ext/Arith.h> // NA_REAL etc. (kept; widely relied on transitively)
#include <array> // simpson weights
#include <cmath> // std::isfinite
#include <concepts> // std::integral, std::floating_point
#include <iterator> // std::iterator_traits
#include <type_traits> // std::is_arithmetic_v
#include <odelia/ode_util.hpp> // to_passive

namespace plant {
namespace util {

struct index {
  index(size_t x_) : x(x_) {}
  size_t check_bounds(size_t size);
  size_t x;
  operator size_t() {return x;}
};

// Inline so the hot per-node competition loop (Species::compute_competition)
// does not pay a cross-TU call for this trivial check (no LTO in this build).
// std::isfinite inlines to a couple of FP instructions; R_FINITE() compiled to
// a cross-library call into libR's R_finite() here, which dominated the K93
// competition loop. Equivalent semantics (NaN/Inf/NA_real_ all -> false) under
// the package's -O2 (no -ffast-math) build.
inline bool is_finite(double x) {
  return std::isfinite(x);
}

// The same test on a scalar carrying derivative layers: the answer picks a
// branch and carries no derivative, so it is taken at the value. Constrained
// away from the built-in numbers only so the overload above keeps them; nothing
// here is a property worth naming.
template <typename T>
  requires (!std::is_arithmetic_v<std::remove_cvref_t<T>>)
bool is_finite(const T& x) {
  return std::isfinite(odelia::util::to_passive(x));
}

void check_length(size_t received, size_t expected);
size_t check_bounds_r(size_t idx, size_t size);

std::vector<double> seq_len(double from, double to, size_t len);

// Use this to be explicit when a deliberate floating point equality test is
// being made (e.g. checking a value is exactly the one we set). Named so the
// intent is clear at the call site; the package's compile flags (-Wall
// -pedantic) do not enable -Wfloat-equal, so no diagnostic suppression is
// needed here.
inline
bool identical(double a, double b) {
  return a == b;
}

// Rcpp converts size_t -> numeric, and I want to be able to preserve
// NA values while doing base 0 to base 1 index conversion (and v.v.).
// This should take the guesswork and remembering out, and should keep
// NA values preserved.
template <class T>
T base_1_to_0(T x) {
  return x - 1;
}

template <class T>
T base_0_to_1(T x) {
  return x + 1;
}

template <class T_from, class T_to>
T_to base_1_to_0(T_from x) {
  return static_cast<T_to>(base_1_to_0<T_from>(x));
}

template <class T_from, class T_to>
T_to base_0_to_1(T_from x) {
  return static_cast<T_to>(base_0_to_1<T_from>(x));
}

// Based on C++11's is_sorted
template <class ForwardIterator>
bool is_sorted(ForwardIterator first, ForwardIterator last) {
  if (first == last)
    return true;

  ForwardIterator next = first;
  while (++next != last) {
    if (*next < *first)
      return false;
    ++first;
  }
  return true;
}
template <class ForwardIterator>
bool is_decreasing(ForwardIterator first, ForwardIterator last) {
  if (first == last)
    return true;

  ForwardIterator next = first;
  while (++next != last) {
    if (*next > *first)
      return false;
    ++first;
  }
  return true;
}

// Wraps Rcpp::stop, which throws and never returns; the attribute lets the
// compiler see that callers (e.g. catch blocks ending in util::stop) do not
// fall through, avoiding spurious -Wreturn-type warnings.
[[noreturn]] void stop(const std::string&);

void warning(const std::string &);

// What the destination iterator holds decides what a serialiser writes: an
// R-facing double takes the rebind value, an active buffer takes the value
// itself. Only double crosses the R boundary, and this is where it converts.
//
// Writes through the iterator rather than returning, because returning meant
// returning by value: at an active scalar that copy registers a tape slot and
// records an operation, once per element per stage, for a value the caller was
// about to assign anyway.
template <typename It, typename S>
void write_iterator_scalar(It& it, const S& value) {
  if constexpr (std::floating_point<
                  typename std::iterator_traits<It>::value_type>) {
    *it++ = odelia::util::to_passive(value);
  } else {
    *it++ = value;
  }
}

// A diagnostic is read rather than differentiated, so an active scalar formats
// at its value. std::to_string covers the built-in types and nothing else.
template<typename T>
std::string to_string(T x) {
  if constexpr (std::integral<T> || std::floating_point<T>) {
    return std::to_string(x);
  } else {
    return std::to_string(odelia::util::to_passive(x));
  }
}

// to_string() always gives six decimal places, so quantities much smaller than
// 1e-6 print as "0.000000" -- useless in an error message about spacings and
// tolerances. This keeps six *significant* figures instead, so 15.7 and 1.5e-05
// both read correctly.
std::string format_double(double x);

template <class ForwardIterator>
void rescale(ForwardIterator first, ForwardIterator last,
             double min_old, double max_old,
             double min_new, double max_new) {
  const double scale = (max_new - min_new) / (max_old - min_old);
  while (first != last) {
    *first = min_new + (*first - min_old) * scale;
    ++first;
  }
}

// Integration via the trapezium rule, for any containers that
// implement the basics of iteration (const_iterator, begin, size)
// The integrand's scalar comes out, so a vector of double integrates to a
// double and nothing about that path changes.
template <typename ContainerX, typename ContainerY>
typename ContainerY::value_type trapezium(const ContainerX& x,
                                          const ContainerY& y) {
  util::check_length(y.size(), x.size());
  if (x.size() < 2) {
    util::stop("Need at least two points for the trapezium rule");
  }
  typename ContainerX::const_iterator x0 = x.begin(), x1 = x.begin();
  ++x1;
  typename ContainerY::const_iterator y0 = y.begin(), y1 = y.begin();
  ++y1;
  typename ContainerY::value_type tot = 0.0;
  while (x1 != x.end()) {
    tot += (*x1++ - *x0++) * (*y1++ + *y0++);
  }
  return tot * 0.5;
}

// 3-point Gauss-Legendre on [-1, 1].
inline constexpr double gauss_legendre_3_x[3] = {-0.7745966692414834, 0.0,
                                                 0.7745966692414834};
inline constexpr double gauss_legendre_3_w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

// Simpson's rule on a non-uniform grid, as three weights on (y0, y1, y2): the
// parabola through (x0, x1, x2) integrated over [x0, x2]. Where either interval
// has zero width the parabola does not exist, and the weights are the
// trapezium's on the two intervals.
inline std::array<double, 3> simpson_pair_weights(double x0, double x1, double x2) {
  const double h0 = x1 - x0, h1 = x2 - x1;
  if (!(h0 > 0.0) || !(h1 > 0.0)) {
    return {h0 / 2, (h0 + h1) / 2, h1 / 2};
  }
  const double s = (h0 + h1) / 6;
  return {s * (2 - h1 / h0), s * (h0 + h1) * (h0 + h1) / (h0 * h1),
          s * (2 - h0 / h1)};
}

// The same parabola integrated over [x1, x2] alone: the last interval of a grid
// with an odd number of intervals.
inline std::array<double, 3> simpson_end_weights(double x0, double x1, double x2) {
  const double h0 = x1 - x0, h1 = x2 - x1;
  if (!(h0 > 0.0) || !(h1 > 0.0)) {
    return {0.0, h1 / 2, h1 / 2};
  }
  return {-h1 * h1 * h1 / (6 * h0 * (h0 + h1)), h1 / 6 * (3 + h1 / h0),
          h1 * (2 * h1 + 3 * h0) / (6 * (h0 + h1))};
}

// Composite Simpson on a non-uniform grid: parabolas over consecutive pairs of
// intervals from the first point, and the last interval of an odd count by
// simpson_end_weights. Exact for quadratics. Two points take the trapezium.
template <typename ContainerX, typename ContainerY>
typename ContainerY::value_type simpson(const ContainerX& x,
                                        const ContainerY& y) {
  util::check_length(y.size(), x.size());
  const size_t n = x.size();
  if (n < 2) {
    util::stop("Need at least two points for Simpson's rule");
  }
  typename ContainerY::value_type tot = 0.0;
  if (n == 2) {
    tot += (x[1] - x[0]) * (y[0] + y[1]) * 0.5;
    return tot;
  }
  size_t i = 0;
  for (; i + 2 < n; i += 2) {
    const std::array<double, 3> w = simpson_pair_weights(x[i], x[i + 1], x[i + 2]);
    tot += w[0] * y[i] + w[1] * y[i + 1] + w[2] * y[i + 2];
  }
  if (i + 1 < n) {
    const std::array<double, 3> w =
      simpson_end_weights(x[n - 3], x[n - 2], x[n - 1]);
    tot += w[0] * y[n - 3] + w[1] * y[n - 2] + w[2] * y[n - 1];
  }
  return tot;
}

template <typename ContainerX, typename ContainerY>
std::vector<double> trapezium_vector(const ContainerX& x,
                                     const ContainerY& y) {
  util::check_length(y.size(), x.size());
  if (x.size() < 2) {
    util::stop("Need at least two points for the trapezium rule");
  }
  typename ContainerX::const_iterator x0 = x.begin(), x1 = x.begin();
  ++x1;
  typename ContainerY::const_iterator y0 = y.begin(), y1 = y.begin();
  ++y1;
  std::vector<double> ret;
  ret.reserve(x.size());
  while (x1 != x.end()) {
    ret.push_back(0.5 * (*x1++ - *x0++) * (*y1++ + *y0++));
  }
  return ret;
}

template <typename T>
T clamp(T x, T min_val, T max_val) {
  return std::max(std::min(x, max_val), min_val);
}

bool is_function(SEXP x);

// The basic idea here is that we consider the three points
//   {(x1, y1), (x2, y2), (x3, y3)}
// and we want to know how much the middle point is contributing to
// the integral, using trapezium integration.
std::vector<double> local_error_integration(const std::vector<double>& x,
                                            const std::vector<double>& y,
                                            double scal);

SEXP get_from_package(const std::string& name);

}
}

namespace Rcpp {
template <> SEXP wrap(const plant::util::index&);
template <> plant::util::index as(SEXP);
template <> SEXP wrap(const std::vector<plant::util::index>&);
}

#endif
