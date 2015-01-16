#include <tree2/util.h>
#include <gsl/gsl_errno.h>
#include <Rcpp.h>
namespace util {

size_t index::check_bounds(size_t size) {
  // We don't check size < 0 or x < 0, as not possible with size_t
  if (size == 0) {
    Rcpp::stop("Index " + util::to_string(x + 1) +
               " out of bounds: container is empty");
  } else if (x >= size) {
    Rcpp::stop("Index " + util::to_string(x + 1) +
               " out of bounds: must be in [1," +
               util::to_string(size) + "]");
  }
  return x;
}

void check_length(size_t received, size_t expected) {
  if (expected != received) {
    Rcpp::stop("Incorrect length input; expected " +
               std::to_string(expected) + ", received " +
               std::to_string(received));
  }
}

bool is_finite(double x) {
  return R_FINITE(x);
}

size_t check_bounds_r(size_t idx, size_t size) {
  // We don't check size < 0 or idx < 0, as not possible with size_t
  if (size == 0) {
    Rcpp::stop("Index " + util::to_string(idx) +
               " out of bounds: container is empty");
  } else if (idx < 1 || idx > size) {
    Rcpp::stop("Index " + util::to_string(idx) +
               " out of bounds: must be in [1," +
               util::to_string(size) + "]");
  }
  return idx - 1;
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

void stop(const std::string& msg) {
  Rcpp::stop(msg);
}

void handler_pass_to_R(const char *reason,
                       const char *file,
                       int line,
                       int gsl_errno) {
  std::ostringstream o;
  o << "GSLERROR: " << reason << ": " <<  file << ":" << line <<
    " [" << util::to_string(gsl_errno) << "]";
  Rcpp::stop(o.str());
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
    for (size_t i = 0; i < x.size(); ++i) {
      ret.push_back(NA_REAL);
    }
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

}

namespace Rcpp {
template <> SEXP wrap(const util::index& x) {
  return Rcpp::wrap(util::base_0_to_1<size_t, int>(x.x));
}
template <> util::index as(SEXP x) {
  const int ix(Rcpp::as<int>(x));
  // TODO: Might need to let NA values through still...
  if (ix <= 0) {
    Rcpp::stop("Invalid value for index (must be >= 1)");
  }
  return util::base_1_to_0<int, size_t>(ix);
}
}

// NOTE: Possibly wants moving?
// [[Rcpp::export]]
void set_sane_gsl_error_handling() {
  gsl_set_error_handler(&util::handler_pass_to_R);
}

// [[Rcpp::export]]
double trapezium(const std::vector<double>& x,
                 const std::vector<double>& y) {
  return util::trapezium(x, y);
}

// [[Rcpp::export]]
std::vector<double> trapezium_vector(const std::vector<double>& x,
                                     const std::vector<double>& y) {
  return util::trapezium_vector(x, y);
}

// [[Rcpp::export]]
std::vector<double> local_error_integration(const std::vector<double>& x,
                                            const std::vector<double>& y,
                                            double scal) {
  return util::local_error_integration(x, y, scal);
}
