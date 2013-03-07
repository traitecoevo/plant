#include <gsl/gsl_errno.h>
#include "util.h"

// See Rcpp-exports for alternative to Rf_error
void handler_pass_to_R(const char *reason,
                       const char *file,
                       int line,
                       int gsl_errno) {
  Rf_error("GSLERROR: %s: %s:%d [%d]", reason, file, line, gsl_errno);
}

// See Rcpp-exports for a better way of doing this.
RcppExport SEXP set_sane_gsl_error_handling() {
  gsl_set_error_handler(&handler_pass_to_R);
  return R_NilValue;
}

