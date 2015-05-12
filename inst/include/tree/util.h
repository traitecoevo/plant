// -*-c++-*-
#ifndef TREE_TREE_UTIL_H_
#define TREE_TREE_UTIL_H_

#include <stddef.h> // size_t
#include <sstream>  // std::stringstream
#include <RcppCommon.h> // as/wrap/SEXP
#include <tree/util_lang_range.h>

namespace tree {
namespace util {

using ::util::lang::indices;
using ::util::lang::range;

struct index {
  index(size_t x_) : x(x_) {}
  size_t check_bounds(size_t size);
  size_t x;
  operator size_t() {return x;}
};

inline std::vector<index> index_vector(const std::vector<size_t> x) {
  std::vector<index> ret;
  ret.reserve(x.size());
  for (size_t i : x) {
    ret.push_back(i);
  }
  return ret;
}

bool is_finite(double x);

void check_length(size_t received, size_t expected);
size_t check_bounds_r(size_t idx, size_t size);
void check_dimensions(size_t recieved_rows, size_t recieved_cols,
                      size_t expected_rows, size_t expected_cols);

std::vector<double> seq_len(double from, double to, size_t len);

// Use this to be explicit when a potentially unsafe floating point
// equality test is being made.  I've disabled the clang warnings
// around this use, while other places warnings will still occur.
inline
bool identical(double a, double b) {
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
#endif
  return a == b;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
}

// http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
bool almost_equal(T x, T y, int ulp) {
  // the machine epsilon has to be scaled to the magnitude of the larger value
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::abs(x - y) <=   std::numeric_limits<T>::epsilon()
    * std::max(std::abs(x), std::abs(y))
    * ulp;
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

void stop(const std::string&);

template<typename T>
std::string to_string(T x) {
  std::ostringstream o;
  if (!(o << x)) {
    stop("String conversion failure");
  }
  return o.str();
}

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

// TODO: Probably move these out to their own file?
// Integration via the trapezium rule, for any containers that
// implement the basics of iteration (const_iterator, begin, size)
template <typename ContainerX, typename ContainerY>
double trapezium(const ContainerX& x, const ContainerY& y) {
  util::check_length(y.size(), x.size());
  if (x.size() < 2) {
    util::stop("Need at least two points for the trapezium rule");
  }
  typename ContainerX::const_iterator x0 = x.begin(), x1 = x.begin();
  ++x1;
  typename ContainerY::const_iterator y0 = y.begin(), y1 = y.begin();
  ++y1;
  double tot = 0.0;
  while (x1 != x.end()) {
    tot += (*x1++ - *x0++) * (*y1++ + *y0++);
  }
  return tot * 0.5;
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
}
}

namespace Rcpp {
template <> SEXP wrap(const tree::util::index&);
template <> tree::util::index as(SEXP);
template <> SEXP wrap(const std::vector<tree::util::index>&);
}

#endif
