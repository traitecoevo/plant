// -*-c++-*-
#ifndef PLANT_PLANT_CANOPY_SHAPE_H_
#define PLANT_PLANT_CANOPY_SHAPE_H_

#include <cmath>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <array>
#include <cstddef>
#include <odelia/ode_util.hpp>

namespace plant {

// How the crown intercepts light. Resolved once per strategy in
// prepare_strategy() (string -> enum), never compared on the hot path.
//
// All except FlatTopBox share the same per-plant competition contribution (the
// smooth Yokozawa leaf-area profile Q); they differ in how a plant's own
// assimilation is computed and how the patch light profile is built:
//   DeepCrown    - assimilation integrated over crown depth against the smooth
//                  light profile: the leaf-area-weighted mean of the (concave)
//                  photosynthetic rate. The original plant behaviour.
//   MeanLight    - integrate the *light* over crown depth to a leaf-area-
//                  weighted mean, then a single photosynthesis evaluation on
//                  that mean light. Partway between DeepCrown and CrownCentre: it
//                  captures the mean light exactly but ignores the curvature of
//                  photosynthesis across the within-crown light distribution.
//   CrownCentre  - identical light profile to DeepCrown, but assimilation is a
//                  single evaluation of the light at the crown centre
//                  (z = H*eta_c) rather than any integral over depth.
//   FlatTopBox   - like CrownCentre for assimilation, but the plant's *competition*
//                  contribution is also collapsed into the thin crown-centre
//                  layer (a hard step: full shade below z = H*eta_c, none above)
//                  instead of the smooth Yokozawa profile. A deliberately naive
//                  variant: it casts shade *incorrectly*, so the patch light
//                  profile is discontinuous and the light-environment spline
//                  cannot be built -- the model does not run. See the vignette.
//   FlatTopSoftBox - a runnable version of FlatTopBox: the step competition is
//                  smoothed into a continuous C1 drop concentrated near the
//                  crown centre (so the light environment can be built), but the
//                  shape is still wrong (box-like, not the gradual Yokozawa
//                  taper). It runs but gives a biased fitness landscape -- the
//                  point being that a *wrong* competition profile yields wrong
//                  evolutionary predictions even when it is numerically fine.
//   PPA          - perfect-plasticity approximation: the patch light profile is
//                  built as a *stepped* function (cumulative leaf area floored
//                  into discrete canopy layers); assimilation then reads that
//                  stepped profile at the crown centre, as CrownCentre does. See
//                  FF16_Environment::compute_environment.
enum class ShadingModel {
  DeepCrown, MeanLight, CrownCentre, FlatTopBox, FlatTopSoftBox, PPA
};

inline ShadingModel shading_model_from_string(const std::string& name) {
  if (name == "deep-crown") {
    return ShadingModel::DeepCrown;
  } else if (name == "mean-light") {
    return ShadingModel::MeanLight;
  } else if (name == "crown-centre") {
    return ShadingModel::CrownCentre;
  } else if (name == "flat-top-box") {
    return ShadingModel::FlatTopBox;
  } else if (name == "flat-top-soft-box") {
    return ShadingModel::FlatTopSoftBox;
  } else if (name == "ppa") {
    return ShadingModel::PPA;
  }
  throw std::invalid_argument("Unknown shading model: " + name);
}

// As above, but an empty string selects the supplied per-strategy default
// (the shared Control default is "" so each strategy picks its own).
inline ShadingModel shading_model_from_string(const std::string& name,
                                              ShadingModel fallback) {
  if (name.empty()) {
    return fallback;
  }
  return shading_model_from_string(name);
}

// Canopy profile used by the FF16/TF24/K93 strategies. The equations follow
// the Yokozawa (1995) foliage-profile model, written in terms of the
// height-normalised coordinate u = z / H:
//
//   q(z, H) = 2 eta (1 - u^eta) u^eta / z
//   Q(z, H) = (1 - u^eta)^2
//   Qp(x, H) = (1 - sqrt(x))^(1 / eta) H
//
// In the solvers these functions are called very frequently from two hot
// paths: crown assimilation quadrature and competition/environment rebuilds.
// The primary q() and Q() methods therefore take the height-normalised ratio
// u = z / H directly, so callers can hoist the z / H division out of inner
// loops. The *_from_height() helpers keep the full-height form available for
// less performance-sensitive code and for reading the original equations.
// initialise() selects an eta-specialised multiplication chain once, and
// caches 1 / eta for Qp(), avoiding repeated generic pow setup where possible
// while preserving the original q(z, H), Q(z, H), and Qp(x, H) semantics. When a
// model can choose eta without changing its intended biology, prefer one of
// the specialised values below (1, 2, 4, 8, 10, 12); other eta values still
// work, but fall back to std::pow().
//
// The profile is templated on the scalar S the coordinate and eta carry.

template <typename S = double>
class CanopyShape {
public:
  using value_type = S;

  CanopyShape()
    : eta_(12.0), eta_inverse_(1.0 / 12.0), eta_c_(eta_c(S(12.0))),
      pow_eta_(&pow_eta_12), leaf_above_(&leaf_above_deep) {
  }

  explicit CanopyShape(S eta) {
    initialise(eta);
  }

  void initialise(S eta, ShadingModel shading_model = ShadingModel::DeepCrown) {
    eta_ = eta;
    eta_inverse_ = 1.0 / eta;
    eta_c_ = eta_c(eta);
    pow_eta_ = select_pow_eta(eta);
    shading_model_ = shading_model;
    // Most models cast shade via the smooth Yokozawa Q (leaf_area_above == Q).
    // FlatTopBox collapses it to a hard step; FlatTopSoftBox to a smoothed step.
    switch (shading_model) {
    case ShadingModel::FlatTopBox:     leaf_above_ = &leaf_above_box;     break;
    case ShadingModel::FlatTopSoftBox: leaf_above_ = &leaf_above_softbox; break;
    default:                           leaf_above_ = &leaf_above_deep;    break;
    }
  }

  // [eqn 11] Fraction of projected leaf area above the height-normalised
  // coordinate u = z / H -- the shading a plant casts at u. Smooth Yokozawa Q
  // for every model except FlatTopBox, which uses a step at the crown centre.
  // Bound once in initialise(), so the competition hot path makes one predicted
  // indirect call with no branch.
  S leaf_area_above(S z_over_height) const {
    return leaf_above_(*this, z_over_height);
  }

  // Undefined at the crown base, where z is 0: use q_from_height, which carries
  // the height the limit there needs.
  S q(S z_over_height, S z) const {
    const S u_eta = pow_eta(z_over_height);
    return 2.0 * eta_ * (1.0 - u_eta) * u_eta / z;
  }

  S q_from_height(S z, S height) const {
    // The 1 / z above is 0 / 0 at the crown base, so take the limit: 0 for every
    // eta above 1, and 2 / height at eta = 1. The light field's lowest knot asks
    // for exactly this, and the crown integral never does.
    if (z <= 0.0) {
      return eta_ == 1.0 ? S(2.0) / height : S(0.0);
    }
    return q(z / height, z);
  }

  // Q(u) and q(z, H) from the single u^eta both need, where u = z / H and
  // height_inverse = 1 / H. q is exactly -dQ/dz, so the second entry is the
  // negative vertical derivative of the first. The Q returned is bit-for-bit the
  // one Q() returns, which is what lets a fused reduction match the value one.
  //
  // FlatTopSoftBox carries its own smoothstep derivative, below. FlatTopBox is a
  // hard step, so its slope is a point mass with no finite value to return -- and
  // that model already cannot build a light environment at all.
  std::pair<S, S> Q_and_q(S z_over_height, S z, S height_inverse) const {
    if (shading_model_ == ShadingModel::FlatTopBox) {
      throw std::runtime_error("The flat-top-box competition profile is a step, "
                               "so it has no vertical slope and no light "
                               "environment can be built from it");
    }
    if (shading_model_ == ShadingModel::FlatTopSoftBox) {
      return Q_and_q_softbox(z_over_height, height_inverse);
    }
    if (z_over_height > 1.0) {
      return {S(0.0), S(0.0)};
    }
    if (z <= 0.0) {
      // The 1 / z below is 0 / 0 at the crown base, so take the limit: 0 for
      // every eta above 1 and 2 / H at eta = 1. The light field's lowest knot
      // asks for exactly this, and the crown integral never does.
      return {Q(z_over_height), eta_ == 1.0 ? 2.0 * height_inverse : S(0.0)};
    }
    const S u_eta = pow_eta(z_over_height);
    const S tmp = 1.0 - u_eta;
    return {tmp * tmp, 2.0 * eta_ * tmp * u_eta / z};
  }


  // Q as a polynomial in w = (z / H)^eta -- for the smooth profile Q = 1 - 2w + w^2
  // -- and w SEPARATES: w = z^eta * H^-eta. So a reduction over crowns at many
  // heights is three running sums over the crowns, one per power of w, and every
  // height reads the same sums. That is what takes a field build from one walk per
  // height to one walk in total.
  //
  // The count is zero for the box profiles, which are piecewise in z / H with an
  // interior breakpoint and so have no such form; FlatTopBox cannot build a field
  // at all. A caller reads the count once per build and takes the walk when it is
  // zero.
  static constexpr std::size_t max_moments = 3;
  std::size_t n_moments() const {
    switch (shading_model_) {
    case ShadingModel::FlatTopBox:
    case ShadingModel::FlatTopSoftBox:
      return 0;
    default:
      return max_moments;
    }
  }

  // The powers of H^-eta one crown contributes. Out-parameters because these are
  // formed once per crown per build and a returned array of active values copies
  // each one, which with a tape active is a recorded operation apiece.
  void crown_moments(const S& height_inverse,
                     std::array<S, max_moments>& out) const {
    out[0] = S(1.0);
    out[1] = pow_eta(height_inverse);
    out[2] = out[1] * out[1];
  }

  // The polynomial's coefficients times the matching powers of z^eta, so that
  // Q(z, H) is this dotted with crown_moments(1 / H).
  //
  // The magnitudes look alarming and are not: z^2eta and H^-2eta are large and
  // small respectively, but a height only ever reads the sums of crowns that reach
  // it, so H >= z wherever the two meet and their product is at most one.
  void height_weights(const S& z, std::array<S, max_moments>& out) const {
    const S t = pow_eta(z);
    out[0] = S(1.0);
    out[1] = -2.0 * t;
    out[2] = t * t;
  }

  // d/dz of the above, which is what the field's slope channel reads from the same
  // sums. z^eta differentiates to eta * z^(eta-1), and z^(eta-1) is 0/0 at the
  // crown base: the limit is 0 for every eta above 1 and 1 at eta = 1, which are
  // the two cases Q_and_q takes by hand.
  void height_weight_slopes(const S& z, std::array<S, max_moments>& out) const {
    const S t = pow_eta(z);
    const S t_over_z = z > 0.0 ? S(t / z) : (eta_ == 1.0 ? S(1.0) : S(0.0));
    out[0] = S(0.0);
    out[1] = -2.0 * eta_ * t_over_z;
    out[2] = 2.0 * eta_ * t * t_over_z;
  }

  S Q(S z_over_height) const {
    if (z_over_height > 1.0) {
      return S(0.0);
    }
    const S tmp = 1.0 - pow_eta(z_over_height);
    return tmp * tmp;
  }

  S Q_from_height(S z, S height) const {
    if (z > height) {
      return S(0.0);
    }
    return Q(z / height);
  }

  S Qp(S x, S height) const {
    return std::pow(1.0 - std::sqrt(x), eta_inverse_) * height;
  }

  // [eqn 12] Crown-centre coordinate u = z / H. Static because the strategies
  // need the same number for their sapwood and conductance terms, and one
  // formula is better than three.
  static S eta_c(S eta) {
    return 1.0 - 2.0 / (1.0 + eta) + 1.0 / (1.0 + 2.0 * eta);
  }

private:
  typedef S (*pow_eta_fn)(S, S);
  typedef S (*leaf_above_fn)(const CanopyShape&, S);

  // u^eta, by THE SAME ALGORITHM AT BOTH SCALARS.
  //
  // ⚠️ Two algorithms here make the DIFFERENTIATED light environment a different
  // quantity from the FORWARD one, and no check can see it because both are
  // plausible. A multiplication chain at double against a correctly-rounded pow
  // at an active scalar disagree in the last bits: at the default eta = 12 the
  // chain is four roundings (u2=u*u; u4=u2*u2; u8=u4*u4; u8*u4), differing on
  // about 82% of u in (0,1) by up to 8 ulp, and Q = (1 - u^eta)^2 cancels near
  // the crown top and turns that into ~4e-9 relative.
  //
  // The derivative that justified the second algorithm is one nothing may read:
  // TF24_Strategy declares `eta` undifferentiable, for exactly the reason the old
  // comment here gave -- u^eta*log(u) is 0*(-inf) at u = 0 and the guard returns a
  // constant zero, "so a row here would be a silently wrong zero". A branch whose
  // only product is a forbidden row is a branch that earns nothing.
  //
  // The guard stays and stays active-only: pow_eta_ may be pow_eta_general, which
  // is pow, and at u = 0 its eta row is that same NaN. At double the chain gives 0
  // there anyway, so leaving the guard off that side keeps the forward numbers
  // bit-identical to what they have always been.
  S pow_eta(S u) const {
    if constexpr (!std::is_same_v<S, double>) {
      if (odelia::util::to_passive(u) <= 0.0) {
        return S(0.0);
      }
    }
    return pow_eta_(u, eta_);
  }

  // Smooth Yokozawa profile -- the correct shading a crown casts.
  static S leaf_above_deep(const CanopyShape& c, S z_over_height) {
    return c.Q(z_over_height);
  }

  // FlatTopBox: all leaf area collapsed into the thin crown-centre layer, so the
  // crown fully shades everything below z = H*eta_c and nothing above. A step.
  static S leaf_above_box(const CanopyShape& c, S z_over_height) {
    return z_over_height < c.eta_c_ ? S(1.0) : S(0.0);
  }

  // FlatTopSoftBox: the hard step softened into a monotone C1 drop, full shade up
  // to lo = max(0, 2*eta_c - 1) then a cubic-smoothstep fall to zero at the crown
  // top (so the transition is centred on the crown centre eta_c and the profile
  // is continuous -- buildable -- but still box-like, not the Yokozawa taper).
  static S leaf_above_softbox(const CanopyShape& c, S z_over_height) {
    const S lo = c.eta_c_ > 0.5 ? 2.0 * c.eta_c_ - 1.0 : S(0.0);
    if (z_over_height <= lo) {
      return S(1.0);
    }
    if (z_over_height >= 1.0) {
      return S(0.0);
    }
    const S t = (z_over_height - lo) / (1.0 - lo);
    return 1.0 - t * t * (3.0 - 2.0 * t);
  }

  // The smoothstep above and its exact negative vertical derivative, from one
  // evaluation of t. d/dt of 1 - t^2(3 - 2t) is -6t(1 - t), so q = -dQ/dz is
  // 6t(1 - t) / ((1 - lo) H), and it vanishes at both ends of the transition.
  std::pair<S, S> Q_and_q_softbox(S z_over_height, S height_inverse) const {
    const S lo = eta_c_ > 0.5 ? 2.0 * eta_c_ - 1.0 : S(0.0);
    if (z_over_height <= lo) {
      return {S(1.0), S(0.0)};
    }
    if (z_over_height >= 1.0) {
      return {S(0.0), S(0.0)};
    }
    const S t = (z_over_height - lo) / (1.0 - lo);
    return {1.0 - t * t * (3.0 - 2.0 * t),
            6.0 * t * (1.0 - t) * height_inverse / (1.0 - lo)};
  }

  static pow_eta_fn select_pow_eta(S eta) {
    if (eta == 1.0) {
      return &pow_eta_1;
    } else if (eta == 2.0) {
      return &pow_eta_2;
    } else if (eta == 4.0) {
      return &pow_eta_4;
    } else if (eta == 8.0) {
      return &pow_eta_8;
    } else if (eta == 10.0) {
      return &pow_eta_10;
    } else if (eta == 12.0) {
      return &pow_eta_12;
    } else {
      return &pow_eta_general;
    }
  }

  static S pow_eta_general(S u, S eta) {
    using std::pow;
    return pow(u, eta);
  }

  static S pow_eta_1(S u, S) {
    return u;
  }

  static S pow_eta_2(S u, S) {
    return u * u;
  }

  static S pow_eta_4(S u, S) {
    const S u2 = u * u;
    return u2 * u2;
  }

  static S pow_eta_8(S u, S) {
    const S u2 = u * u;
    const S u4 = u2 * u2;
    return u4 * u4;
  }

  static S pow_eta_10(S u, S) {
    const S u2 = u * u;
    const S u4 = u2 * u2;
    const S u8 = u4 * u4;
    return u8 * u2;
  }

  static S pow_eta_12(S u, S) {
    const S u2 = u * u;
    const S u4 = u2 * u2;
    const S u8 = u4 * u4;
    return u8 * u4;
  }

  S eta_;
  S eta_inverse_;
  S eta_c_;
  pow_eta_fn pow_eta_;
  leaf_above_fn leaf_above_;
  // Which shape leaf_area_above() is casting. Q_and_q() reads it to refuse the
  // box profiles, whose derivative is not q.
  ShadingModel shading_model_ = ShadingModel::DeepCrown;
};

}

#endif
