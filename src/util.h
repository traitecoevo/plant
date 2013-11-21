// -*-c++-*-
#ifndef TREE_UTIL_H_
#define TREE_UTIL_H_

#include <Rcpp.h>
#include <gsl/gsl_nan.h>
#include <vector>
#include <cstddef>
#include <sstream> // std::stringstream

namespace util {

// This is the same check as in R_ext/Error.h (R 3.0.1) so we should
// agree here.
#if defined(__GNUC__) && __GNUC__ >= 4
#define NORETURN __attribute__((noreturn))
#else
#define NORETURN
#endif

void handler_pass_to_R(const char *reason,
                       const char *file,
                       int line,
                       int gsl_errno) NORETURN;
void set_sane_gsl_error_handling();

bool is_finite(double x);

size_t check_bounds_r(size_t idx, size_t size);
void check_length(size_t received, size_t expected);
void check_dimensions(size_t recieved_rows, size_t recieved_cols,
		      size_t expected_rows, size_t expected_cols);

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

// Wrap up pointers for dealing with passing in and out of R.
// 
// The basic idea is, for some object -- if we are given a raw object
// (not a pointer), we make a a copy of that object and store the
// pointer to it.  When we are destroyed, we'll delete the new memory.
// If we're given a pointer to an object, then we will store that
// pointer, but not clean up memory on exit.
//
// I'm fairly certain that this is implementing a crappy and
// bug-ridden version of something like a smart pointer.  Hopefully
// it's not too bad, and we can swap it back out later if it is.  The
// use case is a bit odd:
//
// Suppose we make an object that *contains* a wrapped pointer.  Then
// on copy of that object:
//
//   1. if the original object is the owner of the information (i.e.,
//      it was given a raw object when initialising the PtrWrapper),
//      the new object will get a copy and be owner too; there will be
//      two fully independent objects.  This logic applies to the all
//      the children too.
// 
//   2. if the original object is not the owner of the information
//      (i.e., it was given a pointer when initialising the
//      PtrWrapper), then the new object only gets a copy of the
//      pointer and won't clean up after itself.  This also applies to
//      all children.
// 
// This probably makes it a little closer to C++11's shared_ptr, but
// because of the structure of how I want to use it (and because of
// the reasons why this matters with passing from R to C via Rcpp
// modules) the reference counting doesn't need anything clever doing
// as we generally have a controlling class.
// 
// This is only meant for use internally within this project.
template <typename T>
class PtrWrapper {
public:
  PtrWrapper(T *obj) : owner(false), ptr(obj)        {}
  PtrWrapper(T  obj) : owner(true),  ptr(new T(obj)) {}
  PtrWrapper(const PtrWrapper &other)
    : owner(other.owner),
      ptr(owner ? new T(*other.ptr) : other.ptr) {}
  PtrWrapper& operator=(PtrWrapper other) {
    using std::swap;
    swap(owner, other.owner);
    swap(ptr,   other.ptr);
    return *this;
  }
  ~PtrWrapper() {
    if (owner)
      delete ptr;
  }
  // See http://stackoverflow.com/a/4421719/1798863
  const T* operator->() const { return ptr; }
  T* operator->() { return ptr; }
  T* get() const { return ptr; }
  bool standalone() const { return owner; }
  // This is only used for testing, and only in Plant.
  // Semantics here: 
  //   - standalones are never equal.
  //   - non-standalones are equal iff they point to the same memory.
  //   - standalones are never equal to non-standalones.
  bool operator==(const PtrWrapper &rhs) const {
    return owner == rhs.owner && (owner || ptr == rhs.ptr);
  }

private:
  bool owner;
  T *ptr;
};

std::vector<double> seq_len(double from, double to, size_t len);

bool identical(double a, double b);
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

std::vector<int> rbinom_multiple(std::vector<int>::iterator it,
				 std::vector<int>::iterator end,
				 double p);

Rcpp::IntegerMatrix to_rcpp_matrix(std::vector< std::vector<int> > x);
std::vector< std::vector<int> > from_rcpp_matrix(Rcpp::IntegerMatrix x);

Rcpp::NumericMatrix to_rcpp_matrix(std::vector< std::vector<double> > x);
std::vector< std::vector<double> > from_rcpp_matrix(Rcpp::NumericMatrix x);

// Integration via the trapezium rule, for any containers that
// implement the basics of iteration (const_iterator, begin, size)
template <typename ContainerX, typename ContainerY>
double trapezium(const ContainerX& x, const ContainerY& y) {
  util::check_length(y.size(), x.size());
  if (x.size() < 2)
    ::Rf_error("Need at least two points for the trapezium rule");
  typename ContainerX::const_iterator x0 = x.begin(), x1 = x.begin();
  x1++;
  typename ContainerY::const_iterator y0 = y.begin(), y1 = y.begin();
  y1++;
  double tot = 0.0;
  while (x1 != x.end())
    tot += (*x1++ - *x0++) * (*y1++ + *y0++);
  return tot * 0.5;
}

template <typename T>
std::string string_from_address(T *x) {
  std::stringstream ss;
  ss << static_cast<const void*>(x);
  return ss.str();
}

// This computes the sum of vectors a and b, so c[i] = a[i] + b[i].
template <typename T>
std::vector<T> sum(const std::vector<T> &a, const std::vector<T> &b) {
  std::vector<T> c(a.size());
  std::transform(a.begin(), a.end(),
		 b.begin(), c.begin(), std::plus<T>());
  return c;
}

std::string rcpp_class_demangle(std::string x);

namespace test {
std::vector<double> test_sum_double(std::vector<double> a,
				    std::vector<double> b);
std::vector<int> test_sum_int(std::vector<int> a,
			      std::vector<int> b);
Rcpp::IntegerMatrix test_to_rcpp_integer_matrix(Rcpp::List x);
Rcpp::List test_from_rcpp_integer_matrix(Rcpp::IntegerMatrix x);

Rcpp::NumericMatrix test_to_rcpp_numeric_matrix(Rcpp::List x);
Rcpp::List test_from_rcpp_numeric_matrix(Rcpp::NumericMatrix x);
}

}

#endif
