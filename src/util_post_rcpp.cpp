#include <tree2/util_post_rcpp.h>

namespace tree2 {
namespace util {

// Given a vector-of-vectors, copy the vector x[i] into the ith
// *column* of an Rcpp matrix.
Rcpp::NumericMatrix to_rcpp_matrix(std::vector< std::vector<double> > x) {
  const size_t n = x.size();
  Rcpp::NumericMatrix ret(static_cast<int>(x.begin()->size()),
                          static_cast<int>(n));
  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < n; ++i)
    it = std::copy(x[i].begin(), x[i].end(), it);
  return ret;
}

}
}
