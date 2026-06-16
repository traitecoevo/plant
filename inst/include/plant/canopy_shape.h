// -*-c++-*-
#ifndef PLANT_PLANT_CANOPY_SHAPE_H_
#define PLANT_PLANT_CANOPY_SHAPE_H_

#include <cmath>

namespace plant {

class CanopyShape {
public:
  CanopyShape()
    : eta_(12.0), pow_eta_(&pow_eta_12) {
  }

  explicit CanopyShape(double eta) {
    initialise(eta);
  }

  void initialise(double eta) {
    eta_ = eta;
    pow_eta_ = select_pow_eta(eta);
  }

  double q(double z, double height) const {
    const double u_eta = pow_eta_(z / height, eta_);
    return 2.0 * eta_ * (1.0 - u_eta) * u_eta / z;
  }

  double Q(double z, double height) const {
    if (z > height) {
      return 0.0;
    }
    const double tmp = 1.0 - pow_eta_(z / height, eta_);
    return tmp * tmp;
  }

  double Qp(double x, double height) const {
    return std::pow(1.0 - std::sqrt(x), 1.0 / eta_) * height;
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
  pow_eta_fn pow_eta_;
};

}

#endif
