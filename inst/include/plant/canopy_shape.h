// -*-c++-*-
#ifndef PLANT_PLANT_CANOPY_SHAPE_H_
#define PLANT_PLANT_CANOPY_SHAPE_H_

#include <cmath>

namespace plant {

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
