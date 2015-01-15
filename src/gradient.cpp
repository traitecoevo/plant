#include <tree2/gradient.h>
#include <tree2/util.h>
#include <tree2/util_post_rcpp.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double test_gradient_fd1(Rcpp::Function f, double x, double dx,
                         int direction, double fx=NA_REAL) {
  if (util::is_finite(fx)) {
    return util::gradient_fd(util::RFunctionWrapper(f), x, dx, fx, direction);
  } else {
    return util::gradient_fd(util::RFunctionWrapper(f), x, dx, direction);
  }
}

// [[Rcpp::export]]
double test_gradient_richardson(Rcpp::Function f, double x, double d,
                                size_t r) {
  return util::gradient_richardson(util::RFunctionWrapper(f), x, d, r);
}
