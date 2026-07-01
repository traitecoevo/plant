// -*-c++-*-
#ifndef PLANT_GRADIENT_COUPLED_CANOPY_H_
#define PLANT_GRADIENT_COUPLED_CANOPY_H_

// Shared coupled-resident canopy kernel (#472 scope B, refactor+optimize phase). The
// Yokozawa light-competition trapezium is identical across FF16 and TF24/TF24f -- in
// all of them area_leaf is pure allometry and competition is density*k_I*area_leaf*Q
// with Q = (1-(z/h)^eta)^2; the leaf optimiser (TF24) affects only the demographic
// RATES, not the light-field geometry. So the
// per-RK-stage canopy reconstruction shared this code by hand-copy (ff16's
// coupled_comp_at == tf24f's tf24f_comp_at); this is the single templated source.

#include <vector>
#include <cstddef>
#include <cmath>
#include <XAD/XAD.hpp>   // active-type value extraction for the crown-height check

namespace plant {
namespace gradient {

// Double value of a scalar that may be a plain double or an XAD active type.
inline double scalar_value(double v) { return v; }
template <typename T> double scalar_value(const T& v) { return xad::value(v); }

// Canopy competition at height z: the descending-height trapezium of geff_i * Q(z/h_i),
// Q = (1-(z/h)^eta)^2, over the active stand (h sorted descending, geff the per-cohort
// weight density*k_I*area_leaf). Cohorts whose crown is below z contribute nothing.
// Returned UN-divided by patch area (the caller applies Beer's law / area). Active in
// the heights/weights, so every trait that moves a height or density re-shades the
// canopy -- the resident coupling. Bit-for-bit the hand-rolled comp_at it replaces.
template <typename S>
S canopy_comp_at(double z, const std::vector<S>& h, const std::vector<S>& geff,
                 double eta) {
  using std::pow;
  const std::size_t n = h.size();
  if (n < 2) return S(0.0);
  auto g = [&](std::size_t i) -> S {
    if (z >= scalar_value(h[i])) return S(0.0);   // no leaf area above the crown
    const S u  = S(z) / h[i];
    const S om = S(1.0) - pow(u, S(eta));
    return geff[i] * (om * om);
  };
  S comp = S(0.0);
  S gp = g(0); S hp = h[0];
  for (std::size_t i = 1; i < n; ++i) {
    S gi = g(i);
    comp = comp + (hp - h[i]) * (gp + gi);
    hp = h[i]; gp = gi;
  }
  return S(0.5) * comp;
}

} // namespace gradient
} // namespace plant

#endif
