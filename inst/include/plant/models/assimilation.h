// -*-c++-*-
#ifndef PLANT_PLANT_ASSIMILATION_H_
#define PLANT_PLANT_ASSIMILATION_H_

#include <memory>
#include <plant/environment.h>
#include <plant/control.h>
#include <plant/qag.h> // quadrature::intervals_type

namespace plant {

template <typename E>
class Assimilation {
  public:

  Assimilation() : a_p1(0.0), a_p2(0.0), eta(0.0) { };
  Assimilation(double a_p1, double a_p2, double eta) : a_p1(a_p1), a_p2(a_p2), eta(eta) { };

  // Leaf productivity parameters
  double a_p1, a_p2;

  // Canopy shape parameter
  double eta;


  // [eqn 12] Gross annual CO2 assimilation
  //
  // NOTE: In contrast with Daniel's implementation (but following
  // Falster 2012), we do not normalise by a_y*a_bio here.
  double assimilate(Control& control,
                    const E& environment,
                    double height,
                    double area_leaf,
                    bool reuse_intervals
                    ) {
    const bool over_distribution = control.plant_assimilation_over_distribution;
    const double x_min = 0.0, x_max = over_distribution ? 1.0 : height;

    double A = 0.0;

    std::function<double(double)> f;
    if (over_distribution) {
      f = [&] (double x) -> double {
        return compute_assimilation_p(x, height, environment);
      };
    } else {
      f = [&] (double x) -> double {
        return compute_assimilation_h(x, height, environment);
      };
    }

    if (control.plant_assimilation_adaptive && reuse_intervals) {
      A = control.integrator.integrate_with_last_intervals(f, x_min, x_max);
    } else {
      A = control.integrator.integrate(f, x_min, x_max);
    }

    return area_leaf * A;
  }

  // This is used in the calculation of assimilation by
  // `compute_assimilation` above; it is the term within the integral in
  // [eqn 12]; i.e., A_lf(A_0v, E(z,a)) * q(z,h(m_l))
  // where `z` is height.
  /* double compute_assimilation_x(Control control, double x, double height, */
  /*                                      const E& environment) const { */
  /*   if (control.plant_assimilation_over_distribution) { */
  /*     return compute_assimilation_p(x, height, environment); */
  /*   } else { */
  /*     return compute_assimilation_h(x, height, environment); */
  /*   } */
  /* } */

  double compute_assimilation_h(double z, double height,
                                       const E& environment) const {
    return assimilation_leaf(environment.get_environment_at_height(z)) * q(z, height);
  }

  double compute_assimilation_p(double p, double height,
                                const E& environment) const {
    return assimilation_leaf(environment.get_environment_at_height(Qp(p, height)));
  }

  // [Appendix S6] Per-leaf photosynthetic rate.
  // Here, `x` is openness, ranging from 0 to 1.
  double assimilation_leaf(double x) const {
    return a_p1 * x / (x + a_p2);
  }

  // (inverse of [eqn 10]; return the height above which fraction 'x' of
  // the leaf mass would be found).
  double Qp(double x, double height) const { // x in [0,1], unchecked.
    return pow(1 - sqrt(x), (1/eta)) * height;
  }

  // [eqn  9] Probability density of leaf area at height `z`
  double q(double z, double height) const {
    const double tmp = pow(z / height, eta);
    return 2 * eta * (1 - tmp) * tmp / z;
  }

  void initialize(double a1, double a2, double e) {
    a_p1 = a1;
    a_p2 = a2;
    eta = e;
  }

};

} // namespace plant

#endif
