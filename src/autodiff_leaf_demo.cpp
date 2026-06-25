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

// --- IFT-through-root-find demo (psi_stem_to_ci) ----------------------------
// Mirrors Leaf::psi_stem_to_ci's residual: colimited assimilation demand minus
// the linear gc supply, solved for ci. The dependence on psi_stem enters via the
// stomatal-conductance coefficient gc, so we differentiate ci w.r.t. gc here
// (dci/dpsi_stem then follows by chain rule with dgc/dpsi_stem, which is the
// spline-based, already-AD-able transport side).
struct CiParams {
  double vcmax = 100.0, et = 120.0, gstar = 4.0, km = 70.0, R_d = 1.0,
         curv = 0.99, ca = 40.0, atm_kPa = 101325.0, umol_to_mol = 1e-6;
};

template <typename T>
T assim_colimited_demo(T ci, const CiParams& p) {
  T ar = p.vcmax * (ci - p.gstar) / (ci + p.km);
  T ae = p.et / 4.0 * (ci - p.gstar) / (ci + 2.0 * p.gstar);
  T s = ar + ae;
  return (s - sqrt(s * s - 4.0 * p.curv * ar * ae)) / (2.0 * p.curv) - p.R_d;
}

// Residual g(ci; gc) whose root defines ci(gc). Templated so AD gives partials.
template <typename T>
T ci_residual(T ci, T gc, const CiParams& p) {
  return assim_colimited_demo(ci, p) * p.umol_to_mol - gc * (p.ca - ci) / p.atm_kPa;
}

// Plain-double bisection solver for ci over (gstar, ca).
double solve_ci(double gc, const CiParams& p) {
  double lo = p.gstar * 1.0001, hi = p.ca * 0.9999;
  for (int it = 0; it < 200; ++it) {
    double mid = 0.5 * (lo + hi);
    if (ci_residual<double>(mid, gc, p) < 0.0) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector ad_ift_demo(double gc = 0.1) {
  CiParams p;
  const double ci_star = solve_ci(gc, p);   // converged root in double

  using AD = xad::fwd<double>::active_type;
  // IFT: dci/dgc = -(dg/dgc)/(dg/dci), partials by forward AD at (ci_star, gc).
  AD ci_a = ci_star; xad::derivative(ci_a) = 1.0;            // seed d/dci
  AD gc_a = gc;
  const double g_ci = xad::derivative(ci_residual(ci_a, gc_a, p));

  AD ci_b = ci_star;
  AD gc_b = gc; xad::derivative(gc_b) = 1.0;                 // seed d/dgc
  const double g_gc = xad::derivative(ci_residual(ci_b, gc_b, p));

  const double dci_dgc_ift = -g_gc / g_ci;

  // Finite-difference check: re-solve the root at gc +/- h.
  const double h = 1e-6 * gc;
  const double dci_dgc_fd = (solve_ci(gc + h, p) - solve_ci(gc - h, p)) / (2 * h);

  return Rcpp::NumericVector::create(
      Rcpp::_["ci_star"] = ci_star,
      Rcpp::_["dci_dgc_ift"] = dci_dgc_ift,
      Rcpp::_["dci_dgc_fd"] = dci_dgc_fd);
}

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
