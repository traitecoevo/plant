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

void set_sane_gsl_error_handling() {
  gsl_set_error_handler(&handler_pass_to_R);
}

size_t check_bounds_r(size_t idx, size_t size) {
  // We don't check size < 0 or idx < 0, as not possible with size_t
  if ( size == 0 )
    ::Rf_error("Index %d out of bounds: container is empty", idx);
  else if ( idx < 1 || idx > size )
    ::Rf_error("Index %d out of bounds: must be in [1,%d]", idx, size);
  return idx - 1;
}

void check_length(size_t received, size_t expected) {
  if ( expected != received )
    ::Rf_error("Incorrect length input; expected %d, received %d\n",
	       expected, received);
}

std::vector<double> seq_len(double from, double to, int len) {
  std::vector<double> ret;
  ret.reserve(len);
  const double dx = (to - from) / (len - 1);
  double x = from;
  for ( int i = 0; i < len; i++, x += dx )
    ret.push_back(x);
  ret.back() = to; // Protect against rounding errors.
  return ret;
}

std::vector<int> rbinom_multiple(std::vector<int>::iterator it,
				 std::vector<int>::iterator end,
				 double p) {
  std::vector<int> ret;
  while ( it != end ) {
    const int k = (int)R::rbinom(*it, p);
    *it -= k;
    ret.push_back(k);
    it++;
  }
  return ret;
}

// Given a vector-of-vectors, copy the vector x[i] into the ith
// *column* of an Rcpp matrix/
Rcpp::IntegerMatrix to_rcpp_matrix(std::vector< std::vector<int> > x) {
  const size_t n = x.size();
  Rcpp::IntegerMatrix ret(x.begin()->size(), n);
  Rcpp::IntegerMatrix::iterator it = ret.begin();
  for ( size_t i = 0; i < n; i++ )
    it = std::copy(x[i].begin(), x[i].end(), it);
  return ret;
}

namespace test {
std::vector<double> test_sum_double(std::vector<double> a,
				    std::vector<double> b) {
  check_length(a.size(), b.size());
  return sum(a, b);
}

std::vector<int> test_sum_int(std::vector<int> a,
			      std::vector<int> b) {
  check_length(a.size(), b.size());
  return sum(a, b);
}

Rcpp::IntegerMatrix test_to_rcpp_matrix(Rcpp::List x) {
  if ( x.size() == 0 )
    Rf_error("Must give positive size 'x'");
  std::vector< std::vector<int> > tmp;
  for ( int i = 0; i < x.size(); i++ ) {
    tmp.push_back(Rcpp::as< std::vector<int> >(x[i]));
  }
  const int n = tmp.begin()->size();
  for ( size_t i = 0; i < tmp.size(); i++ )
    check_length(tmp[i].size(), n);
  return to_rcpp_matrix(tmp);
}
}

}
