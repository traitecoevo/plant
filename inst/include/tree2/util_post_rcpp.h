// -*-c++-*-
#ifndef TREE_UTIL_POST_RCPP_H_
#define TREE_UTIL_POST_RCPP_H_

// Utility bits that need Rcpp.h loaded.  Use this primarily within
// cpp files, not header files.

#include <Rcpp.h>

namespace util {

struct RFunctionWrapper {
  RFunctionWrapper(Rcpp::Function f_) : f(f_) {}
  inline double operator()(double x) {
    return Rcpp::as<double>(f(Rcpp::wrap(x)));
  }
  Rcpp::Function f;
};

}

#endif
