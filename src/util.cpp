#include <tree2/util.h>
#include <Rcpp.h>

namespace util {

void check_length(size_t received, size_t expected) {
  if (expected != received) {
    Rcpp::stop("Incorrect length input; expected " +
               std::to_string(expected) + ", received " +
               std::to_string(received));
  }
}

bool is_finite(double x) {
  return R_FINITE(x);
}

size_t check_bounds_r(size_t idx, size_t size) {
  // We don't check size < 0 or idx < 0, as not possible with size_t
  if (size == 0) {
    Rcpp::stop("Index " + util::to_string(idx) +
               " out of bounds: container is empty");
  } else if (idx < 1 || idx > size) {
    Rcpp::stop("Index " + util::to_string(idx) +
               " out of bounds: must be in [1," +
               util::to_string(size) + "]");
  }
  return idx - 1;
}

}
