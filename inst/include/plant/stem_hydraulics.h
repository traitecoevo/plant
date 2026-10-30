// -*-c++-*-
#ifndef PLANT_PLANT_STEM_HYDRAULICS_H_
#define PLANT_PLANT_STEM_HYDRAULICS_H_

#include <cmath>
#include <odelia/ode_util.hpp>

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
// in TF24_Strategy::net_mass_production_dt).
//
// beta == 0 recovers L_eff = L_top - L_tip, the height-linear model, and with
// L_tip == 0 it returns L_top having performed no arithmetic at all, so the
// collapsed configuration is bit-identical to the old expression rather than
// merely equal to it (invariance criterion I7).
//
// PRECONDITIONS, validated once in TF24_Strategy::prepare_strategy() rather
// than here, because this sits on the compute_rates path:
//   beta == 0  ||  0 < L_tip < L_top
// Templated on the scalar because it sits on the rate path, which the census
// gradient records: L_top is height*eta_c and both L_tip and beta are model
// parameters carrying gradient columns, so all three arrive active on a
// recording pass.
template <typename T>
T effective_path_length(const T& L_top, const T& L_tip, const T& beta) {
  using std::log;
  using std::expm1;
  using odelia::util::to_passive;

  // Both tests read the value and never the derivative. Which closed form the
  // integral takes is piecewise constant in the parameters -- a selector -- and
  // differentiating the choice rather than the model at a fixed choice is what
  // manufactures a discontinuity the model does not have.
  //
  // ⚠️ beta == 0 IS THE ONE CONFIGURATION WHOSE D_c AND theta_c ROWS ARE ZERO
  // RATHER THAN SMALL. The expm1 form below is not merely close to L_top -
  // L_tip at beta == 0, it is that limit exactly, so this arm is a bit-identity
  // shortcut and not a numerical necessity -- but it returns a value that does
  // not read beta, and the recorded row is therefore an exact zero where the
  // true dL_eff/dbeta is L_top*(log(L_top/L_tip) - 1) + L_tip. TF24 defaults to
  // D_c = 0.2, so nothing on the gradient path takes this arm; the height-linear
  // configuration that does (D_c = theta_c = L_tip = 0) exists to assert
  // bit-identity with the pre-v10 model and takes no gradient.
  if (to_passive(beta) == 0.0) {
    // Return the operand verbatim when there is no terminal segment. `L_top -
    // 0.0` is exact in IEEE-754, so this is not about rounding; it removes any
    // dependence on how the compiler treats the subtraction, including whether
    // -ffp-contract fuses the caller's multiply into an fma with it. A
    // regression gate asserting bit-identity should not rest on that reasoning
    // surviving a toolchain upgrade.
    return to_passive(L_tip) == 0.0 ? L_top : T(L_top - L_tip);
  }

  const T x = log(L_top / L_tip);  // > 0 given the precondition
  const T e = T(1.0) - beta;

  // beta == 1 exactly: the logarithmic limit. theta's basipetal decline cancels
  // the conductivity gain term for term, and resistance grows only as log(H).
  if (to_passive(e) == 0.0) {
    return T(L_tip * x);
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
  return T(L_tip * expm1(e * x) / e);
}

}  // namespace stem_hydraulics
}  // namespace plant

#endif
