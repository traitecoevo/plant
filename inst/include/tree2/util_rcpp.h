// -*-c++-*-
#ifndef TREE_UTIL_RCPP_H_
#define TREE_UTIL_RCPP_H_

#include <Rcpp.h>

namespace util {

Rcpp::NumericMatrix to_rcpp_matrix(std::vector< std::vector<double> > x);

}

#endif
