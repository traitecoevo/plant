// -*-c++-*-
#ifndef TREE_TREE_UTIL_POST_RCPP_H_
#define TREE_TREE_UTIL_POST_RCPP_H_

// Utility bits that need Rcpp.h loaded.  Use this primarily within
// cpp files, not header files.

#include <Rcpp.h>

namespace tree {
namespace util {

struct RFunctionWrapper {
  RFunctionWrapper(Rcpp::Function f_) : f(f_) {}
  inline double operator()(double x) {
    return Rcpp::as<double>(f(Rcpp::wrap(x)));
  }
  Rcpp::Function f;
};

Rcpp::NumericMatrix
to_rcpp_matrix(const std::vector< std::vector<double> >& x);
Rcpp::NumericMatrix
to_rcpp_matrix_pad(const std::vector<std::vector<double> >& x);

Rcpp::List matrix_to_list(Rcpp::NumericMatrix x);

}
}

#endif
