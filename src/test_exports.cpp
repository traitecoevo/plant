#include <plant/numerics/fp_contract.h>
#include <plant/models/stem_hydraulics.h>
#include <plant/numerics/gradient.h>
#include <plant/numerics/uniroot.h>
#include <plant/util.h>
#include <plant/util_post_rcpp.h> // RFunctionWrapper
#include <Rcpp.h>

// Test-only entry points, reachable as plant:::test_*() and not
// roxygen-exported, so the unit tests can exercise these functions directly.

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

// The closed-form stem path integral, so that
// tests/testthat/test-tf24-stem-hydraulics.R can exercise the function directly
// without constructing a strategy -- which is what makes the bit-identity
// assertion cheap enough to run over a grid of heights.

// [[Rcpp::export]]
double test_stem_effective_path_length(double L_top, double L_tip, double beta) {
  return plant::stem_hydraulics::effective_path_length(L_top, L_tip, beta);
}

// [[Rcpp::export]]
double test_uniroot(Rcpp::Function f, double min, double max) {
  using namespace plant::util;
  const double tol = 1e-6;
  const size_t max_iterations = 100;
  return uniroot(RFunctionWrapper(f), min, max, tol, max_iterations);
}
