#include <plant/util_post_rcpp.h>

namespace plant {
namespace util {

Rcpp::List matrix_to_list(Rcpp::NumericMatrix x) {
  Rcpp::List ret;
  const size_t nr = static_cast<size_t>(x.nrow());
  for (size_t i = 0; i < nr; ++i) {
    ret.push_back(x(i, Rcpp::_));
  }
  return ret;
}

}
}

// [[Rcpp::export]]
Rcpp::List matrix_to_list(Rcpp::NumericMatrix x) {
  return plant::util::matrix_to_list(x);
}
