#include <plant/gradient.h>
#include <plant/util.h>
#include <plant/util_post_rcpp.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double test_gradient_fd1(Rcpp::Function f, double x, double dx,
                         int direction, double fx=NA_REAL) {
  using namespace plant::util;
  if (is_finite(fx)) {
    return gradient_fd(RFunctionWrapper(f), x, dx, fx, direction);
  } else {
    return gradient_fd(RFunctionWrapper(f), x, dx, direction);
  }
}

// [[Rcpp::export]]
double test_gradient_richardson(Rcpp::Function f, double x, double d,
                                size_t r) {
  using namespace plant::util;
  return gradient_richardson(RFunctionWrapper(f), x, d, r);
}
