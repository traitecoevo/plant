#include <plant/uniroot.h>
#include <plant/util_post_rcpp.h> // RFunctionWrapper
#include <Rcpp.h>

// [[Rcpp::export]]
double test_uniroot(Rcpp::Function f, double min, double max) {
  using namespace plant::util;
  const double tol = 1e-6;
  const int max_iterations = 100;
  return uniroot(RFunctionWrapper(f), min, max, tol, max_iterations);
}
