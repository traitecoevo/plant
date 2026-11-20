#include <plant/numerics/fp_contract.h>
#include <plant/models/stem_hydraulics.h>
#include <plant/numerics/gradient.h>
#include <odelia/root_find.hpp>
#include <plant/util.h>
#include <plant/util_post_rcpp.h> // RFunctionWrapper
#include <Rcpp.h>
#include <plant/models/resource_spline.h>

// Test-only entry points, reachable as plant:::test_*() and not
// roxygen-exported, so the unit tests can exercise these functions directly.

// [[Rcpp::export]]
double test_gradient_fd1(Rcpp::Function f, double x, double dx,
                         int direction, double fx=NA_REAL) {
  using namespace plant::util;
  if (is_finite(fx)) {
    return gradient_fd(RFunctionWrapper(f), x, dx, fx, direction);
  } else {
    return gradient_fd(RFunctionWrapper(f), x, dx, direction);
  }
}

// [[Rcpp::export]]
double test_gradient_richardson(Rcpp::Function f, double x, double d,
                                size_t r) {
  using namespace plant::util;
  return gradient_richardson(RFunctionWrapper(f), x, d, r);
}

// Exposes the closed-form stem path integral, so that
// tests/testthat/test-tf24-stem-hydraulics.R can exercise the function directly
// without constructing a strategy -- which is what makes the bit-identity
// assertion cheap enough to run over a grid of heights.

// [[Rcpp::export]]
double test_stem_effective_path_length(double L_top, double L_tip, double beta) {
  return plant::stem_hydraulics::effective_path_length(L_top, L_tip, beta);
}

// [[Rcpp::export]]
double test_uniroot(Rcpp::Function f, double min, double max) {
  const double tol = 1e-6;
  return odelia::root::toms748(plant::util::RFunctionWrapper(f), min, max,
                               {tol, tol});
}

// The light field's floor, read through the spline every model's field is read
// through: its values at `heights`, at double and at a tangent scalar, and how
// often the floor fired on each path. The tangent field shares the double one's
// differentiated tally, as a rebound environment's field does.
// [[Rcpp::export]]
Rcpp::List test_light_field_floor(std::vector<double> state,
                                  std::vector<double> heights) {
  using tangent = odelia::ode::tangent_scalar<double>;
  plant::ResourceSpline<double> field;
  field.r_init_interpolators(state);
  plant::ResourceSpline<tangent> active;
  active.r_init_interpolators(state);
  active.clamps.differentiated = field.clamps.differentiated;
  std::vector<double> value, value_active;
  for (const double h : heights) {
    value.push_back(field.get_value_at_height(h));
    value_active.push_back(
        odelia::util::to_passive(active.get_value_at_height(tangent(h))));
  }
  const int site = plant::CLAMP_LIGHT_FIELD_FLOOR;
  return Rcpp::List::create(
      Rcpp::_["value"] = value, Rcpp::_["value_active"] = value_active,
      Rcpp::_["forward"] = static_cast<double>(field.clamps.forward[site]),
      Rcpp::_["differentiated"] =
          static_cast<double>((*field.clamps.differentiated)[site]),
      Rcpp::_["forward_at_tangent"] =
          static_cast<double>(active.clamps.forward[site]));
}
