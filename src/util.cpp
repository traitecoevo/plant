#include <gsl/gsl_errno.h>
#include "util.h"

namespace util {

void handler_pass_to_R(const char *reason,
                       const char *file,
                       int line,
                       int gsl_errno) {
  std::ostringstream o;
  o << "GSLERROR: " << reason << ": " <<  file << ":" << line <<
    " [" << util::to_string(gsl_errno) << "]";
  Rcpp::stop(o.str());
}

bool is_finite(double x) {
  return R_FINITE(x);
}

size_t check_bounds_r(size_t idx, size_t size) {
  // We don't check size < 0 or idx < 0, as not possible with size_t
  if (size == 0)
    Rcpp::stop("Index " + util::to_string(idx) +
	       " out of bounds: container is empty");
  else if (idx < 1 || idx > size)
    Rcpp::stop("Index " + util::to_string(idx) +
	       " out of bounds: must be in [1," +
	       util::to_string(size) + "]");
  return idx - 1;
}

void check_length(size_t received, size_t expected) {
  if (expected != received)
    Rcpp::stop("Incorrect length input; expected " +
	       util::to_string(expected) + ", received " +
	       util::to_string(received));
}

void check_dimensions(size_t received_rows, size_t received_cols,
		      size_t expected_rows, size_t expected_cols) {
  if (expected_rows != received_rows)
    Rcpp::stop("Incorrect number of rows; expected " +
	       util::to_string(expected_rows) + ", received " +
	       util::to_string(received_rows));
  if (expected_cols != received_cols)
    Rcpp::stop("Incorrect number of columns; expected " +
	       util::to_string(expected_cols) + ", received " +
	       util::to_string(received_cols));
}

std::vector<double> seq_len(double from, double to, size_t len) {
  std::vector<double> ret;
  ret.reserve(len);
  const double dx = (to - from) / (len - 1);
  double x = from;
  for (size_t i = 0; i < len; ++i, x += dx)
    ret.push_back(x);
  ret.back() = to; // Protect against rounding errors.
  return ret;
}

std::vector<int> rbinom_multiple(std::vector<int>::iterator it,
				 std::vector<int>::iterator end,
				 double p) {
  std::vector<int> ret;
  while (it != end) {
    const int k = static_cast<int>(R::rbinom(*it, p));
    *it -= k;
    ret.push_back(k);
    ++it;
  }
  return ret;
}

// Given a vector-of-vectors, copy the vector x[i] into the ith
// *column* of an Rcpp matrix.
Rcpp::IntegerMatrix to_rcpp_matrix(std::vector< std::vector<int> > x) {
  const size_t n = x.size();
  Rcpp::IntegerMatrix ret(static_cast<int>(x.begin()->size()),
			  static_cast<int>(n));
  Rcpp::IntegerMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < n; ++i)
    it = std::copy(x[i].begin(), x[i].end(), it);
  return ret;
}

// Unpack a matrix column-by-column into a vector of vectors.
std::vector< std::vector<int> > from_rcpp_matrix(Rcpp::IntegerMatrix x) {
  std::vector< std::vector<int> > ret;
  for (int i = 0; i < x.ncol(); ++i) {
    Rcpp::IntegerMatrix::Column x_col_i = x(Rcpp::_, i);
    std::vector<int> xi(x_col_i.begin(), x_col_i.end());
    ret.push_back(xi);
  }
  return ret;
}

// Given a vector-of-vectors, copy the vector x[i] into the ith
// *column* of an Rcpp matrix.
Rcpp::NumericMatrix to_rcpp_matrix(std::vector< std::vector<double> > x) {
  const size_t n = x.size();
  Rcpp::NumericMatrix ret(static_cast<int>(x.begin()->size()),
			  static_cast<int>(n));
  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < n; ++i)
    it = std::copy(x[i].begin(), x[i].end(), it);
  return ret;
}

// Unpack a matrix column-by-column into a vector of vectors.
std::vector< std::vector<double> > from_rcpp_matrix(Rcpp::NumericMatrix x) {
  std::vector< std::vector<double> > ret;
  for (int i = 0; i < x.ncol(); ++i) {
    Rcpp::NumericMatrix::Column x_col_i = x(Rcpp::_, i);
    std::vector<double> xi(x_col_i.begin(), x_col_i.end());
    ret.push_back(xi);
  }
  return ret;
}

// Remove the part saying "Rcpp_" at the beginning so that "Foo" and
// "Rcpp_Foo" end up with the same value.
std::string rcpp_class_demangle(std::string x) {
  size_t pos = x.find("Rcpp_");
  if (pos != std::string::npos) {
    if (pos != 0)
      Rcpp::stop("Cannot determine Rcpp class name");
    else
      x = x.substr(5);
  }
  return x;
}

// The basic idea here is that we consider the three points
//   {(x1, y1), (x2, y2), (x3, y3)}
// and we want to know how much the middle point is contributing to
// the integral.
//
// Because this is normalised against the total error, this could be
// really really badly behaved when this goes towards zero.
std::vector<double> local_error_integration(const std::vector<double>& x,
					    const std::vector<double>& y,
					    double scal) {
  std::vector<double> ret;
  check_length(x.size(), y.size());

  if (x.size() < 3) {
    for (size_t i = 0; i < x.size(); ++i)
      ret.push_back(NA_REAL);
  } else {
    ret.push_back(NA_REAL);
    std::vector<double> a = trapezium_vector(x, y);
    std::vector<double>::const_iterator a1 = a.begin(),
      x1 = x.begin(), y1 = y.begin();
    std::vector<double>::const_iterator a2 = a1+1, x3 = x1+2, y3 = y1+2;
    while (x3 != x.end()) {
      const double a123 = *a1++ + *a2++;
      const double a1_3  = 0.5 * (*y1++ + *y3++) * (*x3++ - *x1++);
      ret.push_back(std::abs(a1_3 - a123) / scal);
    }
    ret.push_back(NA_REAL);
  }
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

Rcpp::IntegerMatrix test_to_rcpp_integer_matrix(Rcpp::List x) {
  if (x.size() == 0)
    Rcpp::stop("Must give positive size 'x'");
  std::vector< std::vector<int> > tmp;
  for (int i = 0; i < x.size(); ++i)
    tmp.push_back(Rcpp::as< std::vector<int> >(x[i]));
  const size_t n = tmp.begin()->size();
  for (size_t i = 0; i < tmp.size(); ++i)
    check_length(tmp[i].size(), n);
  return to_rcpp_matrix(tmp);
}

Rcpp::List test_from_rcpp_integer_matrix(Rcpp::IntegerMatrix x) {
  std::vector< std::vector<int> > res = from_rcpp_matrix(x);
  Rcpp::List ret;
  for (std::vector< std::vector<int> >::iterator it = res.begin();
       it != res.end(); ++it)
    ret.push_back(Rcpp::wrap(*it));
  return ret;
}

Rcpp::NumericMatrix test_to_rcpp_numeric_matrix(Rcpp::List x) {
  if (x.size() == 0)
    Rcpp::stop("Must give positive size 'x'");
  std::vector< std::vector<double> > tmp;
  for (int i = 0; i < x.size(); ++i)
    tmp.push_back(Rcpp::as< std::vector<double> >(x[i]));
  const size_t n = tmp.begin()->size();
  for (size_t i = 0; i < tmp.size(); ++i)
    check_length(tmp[i].size(), n);
  return to_rcpp_matrix(tmp);
}

Rcpp::List test_from_rcpp_numeric_matrix(Rcpp::NumericMatrix x) {
  std::vector< std::vector<double> > res = from_rcpp_matrix(x);
  Rcpp::List ret;
  for (std::vector< std::vector<double> >::iterator it = res.begin();
       it != res.end(); ++it)
    ret.push_back(Rcpp::wrap(*it));
  return ret;
}

}

}
