// -*-c++-*-
#ifndef TREE_QK_H_
#define TREE_QK_H_

#include <vector>
#include <cmath> // std::abs
#include <RcppCommon.h> // SEXP

namespace tree2 {
namespace quadrature {

// Gauss-Kronrod Quadrature -- ported from GSL
class QK {
public:
  QK(size_t rule);
  template <typename Function>
  double integrate(Function f, double a, double b);

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
  std::vector<double> fv1;
  std::vector<double> fv2;
  std::vector<double> xgk, wg, wgk;
};

// Skipping dealing with Functors by just requiring something
// callable using templates.
template <typename Function>
double QK::integrate(Function f, double a, double b) {
  const double center          = 0.5 * (a + b);
  const double half_length     = 0.5 * (b - a);
  const double abs_half_length = std::abs(half_length);
  const double f_center        = f(center);

  double result_gauss = 0;
  double result_kronrod = f_center * wgk[n - 1];
  double result_abs = std::abs(result_kronrod);

  if (n % 2 == 0) {
    result_gauss = f_center * wg[n / 2 - 1];
  }

  for (size_t j = 0; j < (n - 1) / 2; j++) {
    const size_t jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
    const double abscissa = half_length * xgk[jtw];
    const double fval1 = f(center - abscissa);
    const double fval2 = f(center + abscissa);
    const double fsum = fval1 + fval2;
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    result_gauss   += wg[j] * fsum;
    result_kronrod += wgk[jtw] * fsum;
    result_abs     += wgk[jtw] * (std::abs(fval1) + std::abs(fval2));
  }

  for (size_t j = 0; j < n / 2; j++) {
    size_t jtwm1 = j * 2;
    const double abscissa = half_length * xgk[jtwm1];
    const double fval1 = f(center - abscissa);
    const double fval2 = f(center + abscissa);
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

}
}

#endif
