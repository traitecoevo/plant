// PROTOTYPE / DEMO for issue #527 (feature/tf24f-autodiff).
//
// Demonstrates forward-mode automatic differentiation (XAD, via odelia) running
// *inside plant's build* and flowing through (a) the leaf's analytic hydraulic
// cost and (b) a SPLINE, to answer "can splines be made AD-compatible?".
//
// The key trick for the spline: locate the bracketing segment using the plain
// double value of the query (xad::value), then evaluate the segment polynomial
// in the AD scalar so the derivative propagates. Knot data stays double.
//
// Not wired into the model; exposed only as an R-callable demo so we can compare
// the AD derivative against a finite-difference check. Remove before merge.
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <XAD/XAD.hpp>

namespace {

// Templated piecewise-linear spline eval. Segment located by value(u) (a plain
// index lookup), polynomial evaluated in the scalar type U so d/du flows.
template <typename U>
U spline_eval(const std::vector<double>& xs, const std::vector<double>& ys, U u) {
  const double uv = xad::value(u);
  if (uv <= xs.front()) return ys.front() + (u - xs.front()) *
                               ((ys[1] - ys[0]) / (xs[1] - xs[0]));
  const std::size_t n = xs.size();
  if (uv >= xs.back()) return ys[n - 1] + (u - xs[n - 1]) *
                              ((ys[n - 1] - ys[n - 2]) / (xs[n - 1] - xs[n - 2]));
  std::size_t i = 0;
  while (i + 1 < n && xs[i + 1] < uv) ++i;          // locate segment by value
  const double slope = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]);
  return ys[i] + (u - xs[i]) * slope;               // evaluate in AD scalar
}

// Analytic Weibull conductivity proportion, templated (matches Leaf::
// proportion_of_conductivity: exp(-(psi/b)^c)).
template <typename T>
T proportion_of_conductivity(T psi, double b, double c) {
  return exp(-pow(psi / b, c));
}

// Hydraulic cost (matches Leaf::hydraulic_cost_TF) but evaluating the
// conductivity proportion via the SPLINE, so AD must flow through the spline.
template <typename T>
T hydraulic_cost_via_spline(T psi_stem, const std::vector<double>& xs,
                            const std::vector<double>& ys,
                            double g1, double beta2) {
  return g1 * pow(1.0 - spline_eval(xs, ys, psi_stem), beta2);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector ad_leaf_demo(double psi_stem, double b = 3.898245,
                                 double c = 2.680147, double g1 = 7.5,
                                 double beta2 = 1.5) {
  // Build a spline of the analytic conductivity proportion over a psi grid.
  std::vector<double> xs, ys;
  for (double p = 0.0; p <= 12.0 + 1e-9; p += 0.02) {
    xs.push_back(p);
    ys.push_back(std::exp(-std::pow(p / b, c)));
  }

  using mode = xad::fwd<double>;
  using AD = mode::active_type;          // FReal<double> (tapeless forward)

  AD psi = psi_stem;
  xad::derivative(psi) = 1.0;            // seed d/d(psi_stem)
  AD cost = hydraulic_cost_via_spline(psi, xs, ys, g1, beta2);  // through spline

  const double ad_val = xad::value(cost);
  const double ad_der = xad::derivative(cost);

  // Finite-difference check on the same spline-based cost (double path).
  auto cost_d = [&](double p) {
    return hydraulic_cost_via_spline<double>(p, xs, ys, g1, beta2);
  };
  const double h = 1e-6;
  const double fd_der = (cost_d(psi_stem + h) - cost_d(psi_stem - h)) / (2 * h);

  // Also the exact analytic derivative (cost via the analytic Weibull, AD),
  // to show AD-through-spline ~= AD-through-analytic.
  AD psi2 = psi_stem; xad::derivative(psi2) = 1.0;
  AD cost_analytic = g1 * pow(1.0 - proportion_of_conductivity(psi2, b, c), beta2);
  const double ad_analytic_der = xad::derivative(cost_analytic);

  return Rcpp::NumericVector::create(
      Rcpp::_["value"] = ad_val,
      Rcpp::_["ad_deriv_via_spline"] = ad_der,
      Rcpp::_["fd_deriv_via_spline"] = fd_der,
      Rcpp::_["ad_deriv_analytic"] = ad_analytic_der);
}
