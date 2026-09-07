// -*-c++-*-
#ifndef PLANT_PLANT_STEM_HYDRAULICS_H_
#define PLANT_PLANT_STEM_HYDRAULICS_H_

#include <cmath>

namespace plant {
namespace stem_hydraulics {

// Effective hydraulic path length of a stem whose anatomy varies as a power of
// distance-from-apex L. See notes/plan-tf24-height-hydraulics.md sec. 4.
//
// Two within-plant profiles, both anchored at the terminal segment:
//
//   D(L)     = D_tip * (L/L_tip)^D_c            conduit diameter (widening)
//   n_A(L)  propto (D(L)/D_tip)^-2              conduit density (packing limit)
//   k_s(L)   = K_s   * (L/L_tip)^(2*D_c)        sapwood-specific conductivity
//   theta(L) = theta * (L/L_tip)^(-theta_c)     leaf area per sapwood area
//
// The exponent on conductivity is 2*D_c, NOT 4*D_c: the D^4 of Hagen-Poiseuille
// survives only along a single continuous conduit, and at the tissue level a
// conserved lumen fraction means widening is paid for by proportionally fewer
// conduits, halving the effective exponent.
//
// Note the SIGN on theta: theta falls basipetally (so the Huber value 1/theta
// rises), hence -theta_c. A positive theta_c means less leaf area supported per
// unit sapwood as you move down the plant. Getting this backwards inverts the
// compensation the whole exercise is about.
//
// Leaf-specific resistance is R_L = (theta / K_s) * L_eff, where
//
//   L_eff = INT_{L_tip}^{L_top} (L / L_tip)^(-beta) dL,   beta = 2*D_c + theta_c
//
// so a caller wanting a conductance forms K_s * theta / L_eff. L_top is the
// representative flow-path length: for TF24 that is height * eta_c, the
// leaf-area-weighted mean leaf height, not the full height (see the call site
// in src/tf24_strategy.cpp).
//
// beta == 0 recovers L_eff = L_top - L_tip, the height-linear model, and with
// L_tip == 0 it returns L_top having performed no arithmetic at all, so the
// collapsed configuration is bit-identical to the old expression rather than
// merely equal to it (invariance criterion I7).
//
// PRECONDITIONS, validated once in TF24_Strategy::prepare_strategy() rather
// than here, because this sits on the compute_rates path:
//   beta == 0  ||  0 < L_tip < L_top
inline double effective_path_length(double L_top, double L_tip, double beta) {
  if (beta == 0.0) {
    // Return the operand verbatim when there is no terminal segment. `L_top -
    // 0.0` is exact in IEEE-754, so this is not about rounding; it removes any
    // dependence on how the compiler treats the subtraction, including whether
    // -ffp-contract fuses the caller's multiply into an fma with it. A
    // regression gate asserting bit-identity should not rest on that reasoning
    // surviving a toolchain upgrade.
    return L_tip == 0.0 ? L_top : L_top - L_tip;
  }

  const double x = std::log(L_top / L_tip);  // > 0 given the precondition
  const double e = 1.0 - beta;

  // beta == 1 exactly: the logarithmic limit. theta's basipetal decline cancels
  // the conductivity gain term for term, and resistance grows only as log(H).
  if (e == 0.0) {
    return L_tip * x;
  }

  // L_tip * ((L_top/L_tip)^e - 1) / e, written with expm1 rather than as the
  // difference of two powers that sec. 4.4 of the design note prescribes. As
  // e -> 0 both powers approach 1 and differencing them destroys every
  // significant digit (at e = 1e-12 about four survive), whereas expm1(e*x)/e
  // holds full relative precision down to |e| ~ 1e-300. The note's form was
  // chosen to make L_tip = 0 admissible without a division, but that is only
  // ever needed at beta == 0, which returns above.
  //
  // beta > 1 needs no branch: e < 0 makes expm1(e*x) negative too, so the
  // quotient stays positive and tends to L_tip/(beta-1) as L_top -> infinity --
  // the saturating regime of sec. 4.1, where resistance approaches a finite
  // asymptote no matter how tall the plant grows.
  return L_tip * std::expm1(e * x) / e;
}

}  // namespace stem_hydraulics
}  // namespace plant

#endif
