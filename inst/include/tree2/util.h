// -*-c++-*-
#ifndef TREE_UTIL_H_
#define TREE_UTIL_H_

#include <stddef.h> // size_t
#include <sstream>  // std::stringstream
#include <RcppCommon.h>

namespace util {

// Strictly this should be *index* not *count*.
struct count {
  count(size_t x_) : x(x_) {}
  size_t check_bounds(size_t size);
  size_t x;
  operator size_t() {return x;}
};

bool is_finite(double x);

void check_length(size_t received, size_t expected);
size_t check_bounds_r(size_t idx, size_t size);
void check_dimensions(size_t recieved_rows, size_t recieved_cols,
                      size_t expected_rows, size_t expected_cols);

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
// template<class T>
// bool almost_equal(T x, T y, int ulp) {
//   // the machine epsilon has to be scaled to the magnitude of the larger value
//   // and multiplied by the desired precision in ULPs (units in the last place)
//   return std::abs(x - y) <=   std::numeric_limits<T>::epsilon()
//     * std::max(std::abs(x), std::abs(y))
//     * ulp;
// }

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

template<typename T>
std::string to_string(T x) {
  std::ostringstream o;
  if (!(o << x)) {
    // This would be better but required Rcpp.h to be loaded:
    // Rcpp::stop("String conversion failure");
    return "";
  }
  return o.str();
}

}

namespace Rcpp {
template <> SEXP wrap(const util::count&);
template <> util::count as(SEXP);
}

#endif
