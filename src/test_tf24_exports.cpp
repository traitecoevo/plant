#include <plant.h>

// Test-only entry points for TF24's Huber profile, reachable as
// plant:::test_*() and not roxygen-exported. They take tangent derivatives
// through the model's own functions, so the tests can check the two things
// nothing else would catch: the hand-written dmass_sapwood_darea_leaf drifting
// from mass_sapwood, and theta_c's hand-attached column at theta_c = 0.

// Values at double, and the derivative of each by tangent: along the resting
// allometry in leaf area (height following as a_l1 * A^a_l2), and in theta_c at
// fixed leaf area and height.
// [[Rcpp::export]]
Rcpp::List test_tf24_stem_profile(plant::TF24_Model<double> s,
                                  double area_leaf, double sapwood_departure) {
  using tangent = odelia::ode::tangent_scalar<double>;
  using odelia::ode::derivative_along;
  using odelia::ode::set_direction;
  using std::pow;
  s.prepare_strategy();

  const double psi = sapwood_departure;
  const double height = s.pars.a_l1 * pow(area_leaf, s.pars.a_l2);
  const double area_sapwood = s.area_sapwood(area_leaf, psi, height);
  const double area_bark = s.area_bark(area_leaf, height);

  double dmass_sapwood_darea_leaf_tangent;
  {
    const plant::TF24_Model<tangent> t = s.rebind<tangent>();
    tangent A = area_leaf;
    set_direction(A, 1.0);
    const tangent H = t.pars.a_l1 * pow(A, t.pars.a_l2);
    const tangent M = t.mass_sapwood(t.area_sapwood(A, tangent(psi), H), H);
    dmass_sapwood_darea_leaf_tangent = derivative_along(M);
  }

  Rcpp::NumericVector by_theta_c(6);
  {
    plant::TF24_Model<tangent> t = s.rebind<tangent>();
    set_direction(t.pars.theta_c, 1.0);
    const tangent A = area_leaf, H = height, P = psi;
    const tangent as = t.area_sapwood(A, P, H);
    const tangent ab = t.area_bark(A, H);
    by_theta_c[0] = derivative_along(as);
    by_theta_c[1] = derivative_along(t.mass_sapwood(as, H));
    by_theta_c[2] = derivative_along(ab);
    by_theta_c[3] = derivative_along(t.mass_bark(ab, H));
    by_theta_c[4] = derivative_along(t.mass_live_per_area_leaf_at_height(H));
    by_theta_c[5] = derivative_along(t.dmass_sapwood_darea_leaf(A, P));
  }
  by_theta_c.names() = Rcpp::CharacterVector::create(
    "area_sapwood", "mass_sapwood", "area_bark", "mass_bark",
    "mass_live_per_area_leaf_at_height", "dmass_sapwood_darea_leaf");

  return Rcpp::List::create(
    Rcpp::_["height"] = height,
    Rcpp::_["huber_base_factor"] = s.huber_base_factor(height),
    Rcpp::_["area_sapwood"] = area_sapwood,
    Rcpp::_["mass_sapwood"] = s.mass_sapwood(area_sapwood, height),
    Rcpp::_["area_bark"] = area_bark,
    Rcpp::_["mass_bark"] = s.mass_bark(area_bark, height),
    Rcpp::_["mass_live_per_area_leaf_at_height"] =
      s.mass_live_per_area_leaf_at_height(height),
    Rcpp::_["dmass_sapwood_darea_leaf"] =
      s.dmass_sapwood_darea_leaf(area_leaf, psi),
    Rcpp::_["dmass_sapwood_darea_leaf_tangent"] =
      dmass_sapwood_darea_leaf_tangent,
    Rcpp::_["by_theta_c"] = by_theta_c);
}
