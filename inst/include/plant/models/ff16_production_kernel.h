// -*-c++-*-
#ifndef PLANT_PLANT_FF16_PRODUCTION_KERNEL_H_
#define PLANT_PLANT_FF16_PRODUCTION_KERNEL_H_

#include <cmath>   // std::pow; XAD provides pow for active types via ADL

// The elementary pieces of FF16's single-plant net mass production,
// templated on the scalar S so that FF16_Model<S> differentiates
// through them by reverse-mode AD with no special handling. FF16_Model's
// methods delegate to them, at every scalar.

namespace plant {

// Leaf area as a function of height. Carries the
// allometric traits a_l1, a_l2, so it is the entry point through which a trait
// reaches both the mass cascade and the resident competition / light spline.
template <typename S>
S ff16_area_leaf(const S& a_l1, const S& a_l2, const S& height) {
  using std::pow;
  return pow(height / a_l1, 1.0 / a_l2);
}

// Photosynthetic rate per leaf area; x is openness in [0,1]. By const
// reference because it runs at every crown quadrature point, where a copy of an
// active scalar is a tape statement.
template <typename S>
S ff16_assimilation_leaf(const S& a_p1, const S& a_p2, const S& x) {
  return a_p1 * x / (x + a_p2);
}

// Total maintenance respiration (linear in the mass cascade).
// The coarse-root term is appended last so that at the default a_cr1 = 0 (hence
// mass_coarse_root == 0) the sum is bit-identical to the pre-#349 expression.
template <typename S>
S ff16_respiration(const S& mass_leaf, const S& mass_sapwood,
                   const S& mass_bark, const S& mass_root,
                   const S& mass_coarse_root, const S& r_l, const S& r_s,
                   const S& r_b, const S& r_r, const S& r_cr) {
  return r_l * mass_leaf + r_b * mass_bark + r_s * mass_sapwood + r_r * mass_root
    + r_cr * mass_coarse_root;
}

// Total turnover.
template <typename S>
S ff16_turnover(const S& mass_leaf, const S& mass_bark, const S& mass_sapwood,
                const S& mass_root, const S& mass_coarse_root, const S& k_l,
                const S& k_b, const S& k_s, const S& k_r, const S& k_cr) {
  return k_l * mass_leaf + k_b * mass_bark + k_s * mass_sapwood + k_r * mass_root
    + k_cr * mass_coarse_root;
}

// Net production from assimilation/respiration/turnover.
template <typename S>
S ff16_net_production_A(const S& a_bio, const S& a_y, const S& assimilation,
                        const S& respiration, const S& turnover) {
  return a_bio * a_y * (assimilation - respiration) - turnover;
}

}  // namespace plant

#endif
