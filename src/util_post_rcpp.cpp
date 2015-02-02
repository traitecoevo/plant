#include <tree2/util_post_rcpp.h>

namespace tree2 {
namespace util {

// Given a vector-of-vectors, copy the vector x[i] into the ith
// *column* of an Rcpp matrix.
Rcpp::NumericMatrix
to_rcpp_matrix(const std::vector<std::vector<double> >& x) {
  const size_t n = x.size();
  Rcpp::NumericMatrix ret(static_cast<int>(x.begin()->size()),
                          static_cast<int>(n));
  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < n; ++i) {
    it = std::copy(x[i].begin(), x[i].end(), it);
  }
  return ret;
}

// As above, but pads with missing values.
Rcpp::NumericMatrix
to_rcpp_matrix_pad(const std::vector< std::vector<double> >& x) {
  const size_t nc = x.size();
  size_t nr = 0;
  for (auto& i : x) {
    nr = std::max(nr, i.size());
  }
  Rcpp::NumericMatrix ret(static_cast<int>(nr), static_cast<int>(nc));
  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < nc; ++i) {
    const size_t ni = x[i].size();
    it = std::copy_n(x[i].begin(), ni, it);
    it = std::fill_n(it, nr - ni, NA_REAL);
  }
  return ret;
}

}
}
