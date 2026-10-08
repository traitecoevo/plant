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
    diameter_stem(n), mass_root(n), mass_live(n), mass_total(n),
    mass_above_ground(n);

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
    mass_live[i]         = s.mass_live(ml, mb, ms, mr);
    mass_total[i]        = s.mass_total(ml, mb, ms, mh, mr);
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
    Rcpp::_["mass_live"]         = mass_live,
    Rcpp::_["mass_total"]        = mass_total,
    Rcpp::_["mass_above_ground"] = mass_above_ground);
}

} // namespace
} // namespace plant

// [[Rcpp::export]]
Rcpp::List FF16_strategy_expand_allometry(plant::FF16_Strategy s,
                                          Rcpp::NumericVector height,
                                          Rcpp::NumericVector area_heartwood,
                                          Rcpp::NumericVector mass_heartwood) {
  return plant::strategy_expand_allometry_impl(s, height, area_heartwood,
                                               mass_heartwood);
}

// [[Rcpp::export]]
Rcpp::List TF24_strategy_expand_allometry(plant::TF24_Strategy<double> s,
                                          Rcpp::NumericVector height,
                                          Rcpp::NumericVector area_heartwood,
                                          Rcpp::NumericVector mass_heartwood) {
  return plant::strategy_expand_allometry_impl(s, height, area_heartwood,
                                               mass_heartwood);
}

// [[Rcpp::export]]
Rcpp::List TF24f_strategy_expand_allometry(plant::TF24f_Strategy<double> s,
                                          Rcpp::NumericVector height,
                                          Rcpp::NumericVector area_heartwood,
                                          Rcpp::NumericVector mass_heartwood) {
  return plant::strategy_expand_allometry_impl(s, height, area_heartwood,
                                               mass_heartwood);
}
