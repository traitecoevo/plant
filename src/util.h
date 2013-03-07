// -*-c++-*-
#include <Rcpp.h>
#include <gsl/gsl_nan.h>

template <typename T>
bool is_finite(T x) {
  Rf_warning("Requesting finite check, but not yet implemented");
  return true;
}

