// -*-c++-*-
#ifndef TREE_UTIL_H_
#define TREE_UTIL_H_

#include <Rcpp.h>
#include <gsl/gsl_nan.h>
#include <vector>
#include <sstream> // std::stringstream

namespace util {

void set_sane_gsl_error_handling();

template <typename T>
bool is_finite(T x) {
  ::Rf_warning("Requesting finite check; not implemented for this type");
  return true;
}
template <>
bool is_finite(double x);

size_t check_bounds_r(size_t idx, size_t size);
void check_length(size_t received, size_t expected);
void check_dimensions(size_t recieved_rows, size_t recieved_cols,
		      size_t expected_rows, size_t expected_cols);

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
  if ( first == last ) 
    return true;

  ForwardIterator next = first;
  while ( ++next != last ) {
    if ( *next > *first )
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
    if ( owner )
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

std::vector<double> seq_len(double from, double to, int len);

std::vector<int> rbinom_multiple(std::vector<int>::iterator it,
				 std::vector<int>::iterator end,
				 double p);

Rcpp::IntegerMatrix to_rcpp_matrix(std::vector< std::vector<int> > x);
std::vector< std::vector<int> > from_rcpp_matrix(Rcpp::IntegerMatrix x);

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

template <typename T>
std::vector<T> sum(std::vector<T> a, const std::vector<T> b) {
  std::vector<T> c(a.size());
  std::transform(a.begin(), a.end(),
		 b.begin(), c.begin(), std::plus<T>());
  return c;
}

namespace test {
std::vector<double> test_sum_double(std::vector<double> a,
				    std::vector<double> b);
std::vector<int> test_sum_int(std::vector<int> a,
			      std::vector<int> b);
Rcpp::IntegerMatrix test_to_rcpp_matrix(Rcpp::List x);
Rcpp::List test_from_rcpp_matrix(Rcpp::IntegerMatrix x);
}

}

#endif
