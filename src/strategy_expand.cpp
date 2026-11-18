#include <plant.h>

// Allometric expansion of a run's size trajectory, computed with the strategy's
// own C++ allometry functions so the R-side *_expand_state() helpers no longer
// duplicate these formulas (see R/ff16.R, R/tf24.R). The strategy is prepared
// once (resolving eta_c etc.) and the column formulas are evaluated per row;
// the column set and ordering match the historical R expand_state output.

namespace plant {
namespace {

template <typename S>
Rcpp::List strategy_expand_allometry_impl(S s,
                                          const Rcpp::NumericVector& height,
                                          const Rcpp::NumericVector& area_heartwood,
                                          const Rcpp::NumericVector& mass_heartwood) {
  // Resolve derived quantities (eta_c, ...). These are not part of the R-facing
  // strategy list, so a strategy arriving from R is unprepared; prepare here.
  s.prepare_strategy();

  const R_xlen_t n = height.size();
  if (area_heartwood.size() != n || mass_heartwood.size() != n) {
    Rcpp::stop("height, area_heartwood and mass_heartwood must have equal length");
  }
  Rcpp::NumericVector area_leaf(n), mass_leaf(n), area_sapwood(n),
    mass_sapwood(n), area_bark(n), mass_bark(n), area_stem(n),
    diameter_stem(n), mass_root(n), mass_coarse_root(n), mass_live(n),
    mass_total(n), mass_above_ground(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double h  = height[i];
    const double al = s.area_leaf(h);
    const double ml = s.mass_leaf(al);
    const double as = s.area_sapwood(al);
    const double ms = s.mass_sapwood(as, h);
    const double ab = s.area_bark(al);
    const double mb = s.mass_bark(ab, h);
    const double ast = s.area_stem(ab, as, area_heartwood[i]);
    const double mr = s.mass_root(al);
    const double mcr = s.mass_coarse_root(ms);
    const double mh = mass_heartwood[i];

    area_leaf[i]         = al;
    mass_leaf[i]         = ml;
    area_sapwood[i]      = as;
    mass_sapwood[i]      = ms;
    area_bark[i]         = ab;
    mass_bark[i]         = mb;
    area_stem[i]         = ast;
    diameter_stem[i]     = s.diameter_stem(ast);
    mass_root[i]         = mr;
    mass_coarse_root[i]  = mcr;
    mass_live[i]         = s.mass_live(ml, mb, ms, mr, mcr);
    mass_total[i]        = s.mass_total(ml, mb, ms, mh, mr, mcr);
    mass_above_ground[i] = s.mass_above_ground(ml, mb, ms, mh);
  }

  return Rcpp::List::create(
    Rcpp::_["area_leaf"]         = area_leaf,
    Rcpp::_["mass_leaf"]         = mass_leaf,
    Rcpp::_["area_sapwood"]      = area_sapwood,
    Rcpp::_["mass_sapwood"]      = mass_sapwood,
    Rcpp::_["area_bark"]         = area_bark,
    Rcpp::_["mass_bark"]         = mass_bark,
    Rcpp::_["area_stem"]         = area_stem,
    Rcpp::_["diameter_stem"]     = diameter_stem,
    Rcpp::_["mass_root"]         = mass_root,
    Rcpp::_["mass_coarse_root"]  = mass_coarse_root,
    Rcpp::_["mass_live"]         = mass_live,
    Rcpp::_["mass_total"]        = mass_total,
    Rcpp::_["mass_above_ground"] = mass_above_ground);
}

// The crown a strategy's CanopyShape gives, on the unit crown u = z / H: Q, the
// fraction of the leaf area above u (the shade the plant casts, under the
// shading model the strategy resolves from its control), and q, its density in
// u (-dQ/du). Q_and_q at H = 1 reads both from the shape's own profile. eta_c is
// the leaf-area-weighted mean of u.
template <typename S>
Rcpp::List canopy_profile_impl(S s, const Rcpp::NumericVector& u) {
  s.prepare_strategy();
  const R_xlen_t n = u.size();
  Rcpp::NumericVector Q(n), q(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const std::pair<double, double> Qq = s.canopy_shape.Q_and_q(u[i], u[i], 1.0);
    Q[i] = Qq.first;
    q[i] = Qq.second;
  }
  return Rcpp::List::create(
    Rcpp::_["Q"]     = Q,
    Rcpp::_["q"]     = q,
    Rcpp::_["eta_c"] = CanopyShape<double>::eta_c(s.pars.eta));
}

} // namespace
} // namespace plant

// [[Rcpp::export]]
Rcpp::List FF16_canopy_profile(plant::FF16_Model<double> s, Rcpp::NumericVector u) {
  return plant::canopy_profile_impl(s, u);
}

// [[Rcpp::export]]
Rcpp::List TF24_canopy_profile(plant::TF24_Model<double> s, Rcpp::NumericVector u) {
  return plant::canopy_profile_impl(s, u);
}

// [[Rcpp::export]]
Rcpp::List TF24f_canopy_profile(plant::TF24f_Model<double> s, Rcpp::NumericVector u) {
  return plant::canopy_profile_impl(s, u);
}

// [[Rcpp::export]]
Rcpp::List TF24t_canopy_profile(plant::TF24t_Model<double> s, Rcpp::NumericVector u) {
  return plant::canopy_profile_impl(s, u);
}

// [[Rcpp::export]]
Rcpp::List TF24w_canopy_profile(plant::TF24w_Model<double> s, Rcpp::NumericVector u) {
  return plant::canopy_profile_impl(s, u);
}

// [[Rcpp::export]]
Rcpp::List K93_canopy_profile(plant::K93_Model s, Rcpp::NumericVector u) {
  return plant::canopy_profile_impl(s, u);
}

// [[Rcpp::export]]
Rcpp::List ES20_canopy_profile(plant::ES20_Model s, Rcpp::NumericVector u) {
  return plant::canopy_profile_impl(s, u);
}

// K93's size is stem diameter in cm, and its competitive size the basal area in
// cm^2; as FF16's and TF24's stem columns, in m and m^2.
// [[Rcpp::export]]
Rcpp::List K93_strategy_expand_allometry(plant::K93_Model s,
                                         Rcpp::NumericVector height) {
  s.prepare_strategy();
  const R_xlen_t n = height.size();
  Rcpp::NumericVector diameter_stem(n), area_stem(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    diameter_stem[i] = height[i] / 100.0;
    area_stem[i]     = s.size_to_basal_area(height[i]) / 1.0e4;
  }
  return Rcpp::List::create(
    Rcpp::_["area_stem"]     = area_stem,
    Rcpp::_["diameter_stem"] = diameter_stem);
}

// [[Rcpp::export]]
Rcpp::List FF16_strategy_expand_allometry(plant::FF16_Model<double> s,
                                          Rcpp::NumericVector height,
                                          Rcpp::NumericVector area_heartwood,
                                          Rcpp::NumericVector mass_heartwood) {
  return plant::strategy_expand_allometry_impl(s, height, area_heartwood,
                                               mass_heartwood);
}

// [[Rcpp::export]]
Rcpp::List TF24_strategy_expand_allometry(plant::TF24_Model<double> s,
                                          Rcpp::NumericVector height,
                                          Rcpp::NumericVector area_heartwood,
                                          Rcpp::NumericVector mass_heartwood) {
  return plant::strategy_expand_allometry_impl(s, height, area_heartwood,
                                               mass_heartwood);
}

// [[Rcpp::export]]
Rcpp::List TF24f_strategy_expand_allometry(plant::TF24f_Model<double> s,
                                          Rcpp::NumericVector height,
                                          Rcpp::NumericVector area_heartwood,
                                          Rcpp::NumericVector mass_heartwood) {
  return plant::strategy_expand_allometry_impl(s, height, area_heartwood,
                                               mass_heartwood);
}

// [[Rcpp::export]]
Rcpp::List TF24t_strategy_expand_allometry(plant::TF24t_Model<double> s,
                                          Rcpp::NumericVector height,
                                          Rcpp::NumericVector area_heartwood,
                                          Rcpp::NumericVector mass_heartwood) {
  return plant::strategy_expand_allometry_impl(s, height, area_heartwood,
                                               mass_heartwood);
}

// [[Rcpp::export]]
Rcpp::List TF24w_strategy_expand_allometry(plant::TF24w_Model<double> s,
                                          Rcpp::NumericVector height,
                                          Rcpp::NumericVector area_heartwood,
                                          Rcpp::NumericVector mass_heartwood) {
  return plant::strategy_expand_allometry_impl(s, height, area_heartwood,
                                               mass_heartwood);
}
