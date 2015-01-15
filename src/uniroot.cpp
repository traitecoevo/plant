#include <tree2/uniroot.h>
#include <tree2/util_post_rcpp.h> // RFunctionWrapper
#include <Rcpp.h>

// [[Rcpp::export]]
double test_uniroot(Rcpp::Function f, double min, double max) {
  const double tol = 1e-6;
  const int max_iterations = 100;
  return util::uniroot(util::RFunctionWrapper(f), min, max, tol,
                       max_iterations);
}
