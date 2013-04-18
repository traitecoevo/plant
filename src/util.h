// -*-c++-*-
#ifndef TREE_UTIL_H_
#define TREE_UTIL_H_

#include <Rcpp.h>
#include <gsl/gsl_nan.h>
#include <vector>
#include <sstream> // std::stringstream

namespace util {

template <typename T>
bool is_finite(T x) {
  // TODO: Get the finite check in here!
  // Rf_warning("Requesting finite check, but not yet implemented");
  return true;
}

void check_bounds(size_t idx, size_t max);
void check_length(size_t received, size_t expected);

// Based on C++11's is_sorted
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

// I'm fairly certain that this is implementing a crappy and
// bug-ridden version of something like a smart pointer.  Hopefully
// it's not too bad, and we can swap it back out later.
// 
// This is only meant for use internally within this project, so
// hopfully nothing else copies it too much.
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
  bool standalone() const { return owner; }
  // Semantics here: 
  //   - standalones are never equal.
  //   - non-standalones are equal iff they point to the same memory.
  //   - standalones are never equal to non-standalones.
  bool operator==(const PtrWrapper &rhs) const {
    return owner == rhs.owner && (owner || ptr == rhs.ptr);
  }

  // public for now...
  bool owner;
  T *ptr;
};

std::vector<double> seq_len(double from, double to, int len);

template <typename T>
std::string string_from_address(T *x) {
  std::stringstream ss;
  ss << static_cast<const void*>(x);
  return ss.str();
}

}

#endif
