#include <plant/qk.h>
#include <plant/qk_rules.h>

#include <plant/util.h> // check_length
#include <plant/util_post_rcpp.h> // RFunctionWrapper

namespace plant {
namespace quadrature {

QK::QK()
  : last_result(NA_REAL),
    last_result_abs(NA_REAL),
    last_result_asc(NA_REAL),
    last_error(NA_REAL) {
  initialise(21);
}


QK::QK(size_t rule)
  : last_result(NA_REAL),
    last_result_abs(NA_REAL),
    last_result_asc(NA_REAL),
    last_error(NA_REAL) {
  initialise(rule);
}

// X points used by the quadrature routine.
std::vector<double> QK::integrate_vector_x(double a, double b) const {
  const double center          = 0.5 * (a + b);
  const double half_length     = 0.5 * (b - a);

  std::vector<double> x;
  x.push_back(center);

  for (size_t j = 0; j < (n - 1) / 2; j++) {
    const size_t jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
    const double abscissa = half_length * xgk[jtw];
    x.push_back(center - abscissa);
    x.push_back(center + abscissa);
  }

  for (size_t j = 0; j < n / 2; j++) {
    size_t jtwm1 = j * 2;
    const double abscissa = half_length * xgk[jtwm1];
    x.push_back(center - abscissa);
    x.push_back(center + abscissa);
  }

  return x;
}

double QK::integrate_vector(const std::vector<double>& y,
                            double a, double b) {
  util::check_length(y.size(), 2 * n - 1);
  std::vector<double>::const_iterator yy = y.begin();

  const double half_length     = 0.5 * (b - a);
  const double abs_half_length = std::abs(half_length);
  const double f_center        = *yy++;

  double result_gauss = 0;
  double result_kronrod = f_center * wgk[n - 1];
  double result_abs = std::abs(result_kronrod);

  if (n % 2 == 0) {
    result_gauss = f_center * wg[n / 2 - 1];
  }

  for (size_t j = 0; j < (n - 1) / 2; j++) {
    const size_t jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
    const double fval1 = *yy++;
    const double fval2 = *yy++;
    const double fsum = fval1 + fval2;
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    result_gauss   += wg[j] * fsum;
    result_kronrod += wgk[jtw] * fsum;
    result_abs     += wgk[jtw] * (std::abs(fval1) + std::abs(fval2));
  }

  for (size_t j = 0; j < n / 2; j++) {
    size_t jtwm1 = j * 2;
    const double fval1 = *yy++;
    const double fval2 = *yy++;
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    result_kronrod += wgk[jtwm1] * (fval1 + fval2);
    result_abs     += wgk[jtwm1] * (std::abs(fval1) + std::abs(fval2));
  };

  const double mean = result_kronrod * 0.5;

  double result_asc = wgk[n - 1] * std::abs(f_center - mean);

  for (size_t j = 0; j < n - 1; j++) {
    result_asc += wgk[j] * (std::abs(fv1[j] - mean) +
                            std::abs(fv2[j] - mean));
  }

  /* scale by the width of the integration region */

  const double err = (result_kronrod - result_gauss) * half_length;

  result_kronrod *= half_length;
  result_abs *= abs_half_length;
  result_asc *= abs_half_length;

  last_result     = result_kronrod;
  last_result_abs = result_abs;
  last_result_asc = result_asc;
  last_error      = rescale_error(err, result_abs, result_asc);

  return last_result;
}

double QK::rescale_error(double err, double result_abs,
                         double result_asc) {
  err = std::abs(err);

  if (result_asc > 0 && err > 0) {
    double scale = pow((200 * err / result_asc), 1.5);
    if (scale < 1) {
      err = result_asc * scale;
    } else {
      err = result_asc;
    }
  }

  const double e = std::numeric_limits<double>::epsilon();
  const double u = std::numeric_limits<double>::min();
  if (result_abs > u / (50 * e)) {
    err = std::max(50 * e * result_abs, err);
  }

  return err;
}

double QK::r_integrate(SEXP f, double a, double b) {
  util::RFunctionWrapper fw(Rcpp::as<Rcpp::Function>(f));
  return integrate(fw, a, b);
}

// This could probably be done way better by providing a base class
// from which of the different rules inherits, but that just seems
// overkill.
void QK::initialise(size_t rule) {
  if (rule == 15) {
    n = QK15::n;
    const size_t ng = QK15::ng;
    xgk.resize(n);
    wg.resize(ng);
    wgk.resize(n);
    std::copy(QK15::xgk, QK15::xgk + n,  xgk.begin());
    std::copy(QK15::wg,  QK15::wg  + ng,  wg.begin());
    std::copy(QK15::wgk, QK15::wgk + n,  wgk.begin());
  } else if (rule == 21) {
    n = QK21::n;
    const size_t ng = QK21::ng;
    xgk.resize(n);
    wg.resize(ng);
    wgk.resize(n);
    std::copy(QK21::xgk, QK21::xgk + n,  xgk.begin());
    std::copy(QK21::wg,  QK21::wg  + ng,  wg.begin());
    std::copy(QK21::wgk, QK21::wgk + n,  wgk.begin());
  } else if (rule == 31) {
    n = QK31::n;
    const size_t ng = QK31::ng;
    xgk.resize(n);
    wg.resize(ng);
    wgk.resize(n);
    std::copy(QK31::xgk, QK31::xgk + n,  xgk.begin());
    std::copy(QK31::wg,  QK31::wg  + ng,  wg.begin());
    std::copy(QK31::wgk, QK31::wgk + n,  wgk.begin());
  } else if (rule == 41) {
    n = QK41::n;
    const size_t ng = QK41::ng;
    xgk.resize(n);
    wg.resize(ng);
    wgk.resize(n);
    std::copy(QK41::xgk, QK41::xgk + n,  xgk.begin());
    std::copy(QK41::wg,  QK41::wg  + ng,  wg.begin());
    std::copy(QK41::wgk, QK41::wgk + n,  wgk.begin());
  } else if (rule == 51) {
    n = QK51::n;
    const size_t ng = QK51::ng;
    xgk.resize(n);
    wg.resize(ng);
    wgk.resize(n);
    std::copy(QK51::xgk, QK51::xgk + n,  xgk.begin());
    std::copy(QK51::wg,  QK51::wg  + ng,  wg.begin());
    std::copy(QK51::wgk, QK51::wgk + n,  wgk.begin());
  } else if (rule == 61) {
    n = QK61::n;
    const size_t ng = QK61::ng;
    xgk.resize(n);
    wg.resize(ng);
    wgk.resize(n);
    std::copy(QK61::xgk, QK61::xgk + n,  xgk.begin());
    std::copy(QK61::wg,  QK61::wg  + ng,  wg.begin());
    std::copy(QK61::wgk, QK61::wgk + n,  wgk.begin());
  } else {
    Rcpp::stop("Unknown rule " + util::to_string(rule));
  }
  fv1.resize(n);
  fv2.resize(n);
}

}
}
