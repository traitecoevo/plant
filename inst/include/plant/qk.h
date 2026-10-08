// -*-c++-*-

// The class `QK` provides methods to use Gauss-Kronrod quadrature to numerically integrate a function.
// This is code is ported from GSL
// The integration has several "rules", defined in qk_rules.cpp. These include QK15, QK21, QK31, QK41, QK51. These allow for different numbers of points in the integration.
// General background on the method is available at: https://en.wikipedia.org/wiki/ or
// or https://www.gnu.org/software/gsl/doc/html/integration.html
// These methods can also be accessed via the class QAG, which offers an adaptive extension of the method.
// Implemented by Rich FitzJohn

#ifndef PLANT_PLANT_QK_H_
#define PLANT_PLANT_QK_H_

#include <vector>
#include <cmath> // std::abs, std::isnan
#include <RcppCommon.h> // SEXP
#include <odelia/ode_util.hpp> // to_passive
#include <plant/util.h> // util::stop

namespace plant {
namespace quadrature {

class QK {
public:
  QK();
  QK(size_t rule);
  // The limits and the integrand's values carry the scalar the caller evaluates
  // at; the rule's abscissae and weights are fractions of the interval and stay
  // double. A crown integral's upper limit is the cohort's height, so reading
  // the field at the fraction instead of at the position built from it makes
  // d(area)/d(height) exactly zero.
  template <typename Function, typename S>
  S integrate(Function f, const S& a, const S& b);

  // These two provide very low level access to the integration
  // routines.
  std::vector<double> integrate_vector_x(double a, double b) const;
  double integrate_vector(const std::vector<double>& y,
                          double a, double b);

  double get_last_area()     const {return last_result;    }
  double get_last_error()    const {return last_error;     }
  double get_last_area_abs() const {return last_result_abs;}
  double get_last_area_asc() const {return last_result_asc;}

  // * R interface
  double r_integrate(SEXP f, double a, double b);

private:
  void initialise(size_t rule);
  static double rescale_error(double err, double result_abs,
                              double result_asc);

  double last_result;
  double last_result_abs;
  double last_result_asc;
  double last_error;

  size_t n;
  // Integrand values at the abscissae, kept only to form last_result_asc, so
  // they hold the values even when the integrand carries a scalar.
  std::vector<double> fv1;
  std::vector<double> fv2;
  std::vector<double> xgk, wg, wgk;
};

// Skipping dealing with Functors by just requiring something
// callable using templates.
template <typename Function, typename S>
S QK::integrate(Function f, const S& a, const S& b) {
  using odelia::util::to_passive;
  const S center               = 0.5 * (a + b);
  const S half_length          = 0.5 * (b - a);
  const double abs_half_length = std::abs(to_passive(half_length));
  const S f_center             = f(center);
  if (std::isnan(to_passive(f_center))) {
    util::stop("Integrand returned NaN at x=" +
               std::to_string(to_passive(center)));
  }

  // The Gauss sum exists only to form the error estimate, which is a last_*
  // diagnostic -- so it is accumulated at the values, like the other three. DO NOT
  // make it S: at the active scalar it records a statement per Gauss node whose
  // row nothing ever reads, because the only thing downstream of it is
  // to_passive(err). That was 8-11% of this rule's tape.
  double result_gauss = 0.0;
  S result_kronrod = f_center * wgk[n - 1];
  // The four last_* diagnostics are read by a user, not differentiated, so the
  // absolute and ascending sums are accumulated at the values throughout.
  double result_abs = std::abs(to_passive(result_kronrod));

  if (n % 2 == 0) {
    result_gauss = to_passive(f_center) * wg[n / 2 - 1];
  }

  for (size_t j = 0; j < (n - 1) / 2; j++) {
    const size_t jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
    const S abscissa = half_length * xgk[jtw];
    const S fval1 = f(center - abscissa);
    const S fval2 = f(center + abscissa);
    const S fsum = fval1 + fval2;
    fv1[jtw] = to_passive(fval1);
    fv2[jtw] = to_passive(fval2);
    result_gauss   += wg[j] * to_passive(fsum);
    result_kronrod += wgk[jtw] * fsum;
    result_abs     += wgk[jtw] * (std::abs(fv1[jtw]) + std::abs(fv2[jtw]));
  }

  for (size_t j = 0; j < n / 2; j++) {
    size_t jtwm1 = j * 2;
    const S abscissa = half_length * xgk[jtwm1];
    const S fval1 = f(center - abscissa);
    const S fval2 = f(center + abscissa);
    fv1[jtwm1] = to_passive(fval1);
    fv2[jtwm1] = to_passive(fval2);
    result_kronrod += wgk[jtwm1] * (fval1 + fval2);
    result_abs     += wgk[jtwm1] * (std::abs(fv1[jtwm1]) + std::abs(fv2[jtwm1]));
  };

  // The mean only centres the ascending sum, which is a diagnostic, so it is
  // taken at the value for the same reason the Gauss sum is.
  const double mean_value = to_passive(result_kronrod) * 0.5;

  double result_asc =
    wgk[n - 1] * std::abs(to_passive(f_center) - mean_value);

  for (size_t j = 0; j < n - 1; j++) {
    result_asc += wgk[j] * (std::abs(fv1[j] - mean_value) +
                            std::abs(fv2[j] - mean_value));
  }

  /* scale by the width of the integration region */

  const double err =
    (to_passive(result_kronrod) - result_gauss) * to_passive(half_length);

  result_kronrod *= half_length;
  result_abs *= abs_half_length;
  result_asc *= abs_half_length;

  last_result     = to_passive(result_kronrod);
  last_result_abs = result_abs;
  last_result_asc = result_asc;
  last_error      = rescale_error(err, result_abs, result_asc);

  if (std::isnan(last_result)) {
    util::stop("Integrand produced NaN result over [" +
               std::to_string(to_passive(a)) + ", " +
               std::to_string(to_passive(b)) + "]");
  }

  return result_kronrod;
}

}
}

#endif
