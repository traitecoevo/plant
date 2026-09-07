#include <plant/stem_hydraulics.h>
#include <Rcpp.h>

// Test-only entry point for the closed-form stem path integral. Reachable as
// plant:::test_stem_effective_path_length(); `test_`-prefixed and not
// roxygen-exported, following test_gradient_fd1 in src/gradient.cpp, so the
// unit tests in tests/testthat/test-tf24-stem-hydraulics.R can exercise the
// function directly without constructing a strategy -- which is what makes the
// bit-identity assertion cheap enough to run over a grid of heights.

// [[Rcpp::export]]
double test_stem_effective_path_length(double L_top, double L_tip, double beta) {
  return plant::stem_hydraulics::effective_path_length(L_top, L_tip, beta);
}
