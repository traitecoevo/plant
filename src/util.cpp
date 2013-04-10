#include <gsl/gsl_errno.h>
#include "util.h"

namespace util {

// See Rcpp-exports for alternative to Rf_error
void handler_pass_to_R(const char *reason,
                       const char *file,
                       int line,
                       int gsl_errno) {
  ::Rf_error("GSLERROR: %s: %s:%d [%d]", reason, file, line, gsl_errno);
}

// See Rcpp-exports for a better way of doing this.
// Actually, do this via Rcpp::function in interface.
RcppExport SEXP set_sane_gsl_error_handling() {
  gsl_set_error_handler(&handler_pass_to_R);
  return R_NilValue;
}

void check_bounds(int idx, int max) {
  if ( max < 0 )
    ::Rf_error("Impossible upper bound");
  else if ( max == 0 )
    ::Rf_error("Index %d impossible from empty range", idx);
  else if ( idx < 0 || idx >= max )
    ::Rf_error("Index %d out of bounds: must be in [0,%d]", idx, max-1);
}

void check_length(size_t received, size_t expected) {
  if ( expected != received )
    ::Rf_error("Incorrect length input; expected %d, received %d\n",
	       expected, received);
}

}
