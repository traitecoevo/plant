// -*-c++-*-
#ifndef PLANT_PLANT_CANOPY_SHAPE_H_
#define PLANT_PLANT_CANOPY_SHAPE_H_

#include <cmath>
#include <string>
#include <stdexcept>

namespace plant {

// How the crown intercepts light. Resolved once per strategy in
// prepare_strategy() (string -> enum), never compared on the hot path.
//
// All share the same per-plant competition contribution (the smooth Yokozawa
// leaf-area profile Q); they differ in how a plant's own assimilation is
// computed and how the patch light profile is built:
//   DeepCrown    - assimilation integrated over crown depth against the smooth
//                  light profile: the leaf-area-weighted mean of the (concave)
//                  photosynthetic rate. The original plant behaviour.
//   AverageLight - integrate the *light* over crown depth to a leaf-area-
//                  weighted mean, then a single photosynthesis evaluation on
//                  that mean light. Partway between DeepCrown and FlatTop: it
//                  captures the mean light exactly but ignores the curvature of
//                  photosynthesis across the within-crown light distribution.
//   FlatTop      - identical light profile to DeepCrown, but assimilation is a
//                  single evaluation of the light at the crown centre
//                  (z = H*eta_c) rather than any integral over depth.
//   PPA          - perfect-plasticity approximation: the patch light profile is
//                  built as a *stepped* function (cumulative leaf area floored
//                  into discrete canopy layers); assimilation then reads that
//                  stepped profile at the crown centre, as FlatTop does. See
//                  FF16_Environment::compute_environment.
enum class ShadingModel { DeepCrown, AverageLight, FlatTop, PPA };

inline ShadingModel shading_model_from_string(const std::string& name) {
  if (name == "deep-crown") {
    return ShadingModel::DeepCrown;
  } else if (name == "average-light") {
    return ShadingModel::AverageLight;
  } else if (name == "flat-top") {
    return ShadingModel::FlatTop;
  } else if (name == "ppa") {
    return ShadingModel::PPA;
  }
  throw std::invalid_argument("Unknown shading model: " + name);
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

class CanopyShape {
public:
  CanopyShape()
    : eta_(12.0), eta_inverse_(1.0 / 12.0), pow_eta_(&pow_eta_12) {
  }

  explicit CanopyShape(double eta) {
    initialise(eta);
  }

  void initialise(double eta) {
    eta_ = eta;
    eta_inverse_ = 1.0 / eta;
    pow_eta_ = select_pow_eta(eta);
  }

  double q(double z_over_height, double z) const {
    const double u_eta = pow_eta_(z_over_height, eta_);
    return 2.0 * eta_ * (1.0 - u_eta) * u_eta / z;
  }

  double q_from_height(double z, double height) const {
    return q(z / height, z);
  }

  double Q(double z_over_height) const {
    if (z_over_height > 1.0) {
      return 0.0;
    }
    const double tmp = 1.0 - pow_eta_(z_over_height, eta_);
    return tmp * tmp;
  }

  double Q_from_height(double z, double height) const {
    if (z > height) {
      return 0.0;
    }
    return Q(z / height);
  }

  double Qp(double x, double height) const {
    return std::pow(1.0 - std::sqrt(x), eta_inverse_) * height;
  }

private:
  typedef double (*pow_eta_fn)(double, double);

  static pow_eta_fn select_pow_eta(double eta) {
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

  static double pow_eta_general(double u, double eta) {
    return std::pow(u, eta);
  }

  static double pow_eta_1(double u, double) {
    return u;
  }

  static double pow_eta_2(double u, double) {
    return u * u;
  }

  static double pow_eta_4(double u, double) {
    const double u2 = u * u;
    return u2 * u2;
  }

  static double pow_eta_8(double u, double) {
    const double u2 = u * u;
    const double u4 = u2 * u2;
    return u4 * u4;
  }

  static double pow_eta_10(double u, double) {
    const double u2 = u * u;
    const double u4 = u2 * u2;
    const double u8 = u4 * u4;
    return u8 * u2;
  }

  static double pow_eta_12(double u, double) {
    const double u2 = u * u;
    const double u4 = u2 * u2;
    const double u8 = u4 * u4;
    return u8 * u4;
  }

  double eta_;
  double eta_inverse_;
  pow_eta_fn pow_eta_;
};

}

#endif
