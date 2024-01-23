// -*-c++-*-
#ifndef PLANT_PLANT_ASSIMILATION_H_
#define PLANT_PLANT_ASSIMILATION_H_

#include <memory>
#include <plant/environment.h>
#include <plant/control.h>
#include <plant/qag.h>

// Used for integrating assimilation over a range of depths within the plant

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

  //Gross annual CO2 assimilation
  double assimilate(const E& environment,
                    double height,
                    double area_leaf,
                    bool reuse_intervals) {

    const double x_min = 0.0, x_max = height;

    double A = 0.0;

    std::function<double(double)> f_assimilation_at_height = [&](double x) -> double
    {
      return  assimilation_leaf(environment.get_environment_at_height(x)) * q(x, height);
    };

    if (integrator.is_adaptive() && reuse_intervals) {
      A = integrator.integrate_with_last_intervals(f_assimilation_at_height, x_min, x_max);
    } else {
      A = integrator.integrate(f_assimilation_at_height, x_min, x_max);
    }

    return area_leaf * A;
  }

  // [Appendix S6] Per-leaf photosynthetic rate.
  // Here, `x` is openness, ranging from 0 to 1.
  double assimilation_leaf(double x) const {
    return a_p1 * x / (x + a_p2);
  }

  // [eqn  9] Probability density of leaf area at height `z`
  double q(double z, double height) const {
    const double tmp = pow(z / height, eta);
    return 2 * eta * (1 - tmp) * tmp / z;
  }

  void initialize(double a1, double a2, double e,
                  bool adaptive_integration=true,
                  int integration_rule=21,
                  int iterations=1000,
                  double integration_tol=1e-6) {

    // strategy parameters
    a_p1 = a1;
    a_p2 = a2;
    eta = e;

    // set up integrator
    if(!adaptive_integration) {
      iterations = 1;
    }

    integrator = quadrature::QAG(integration_rule,
                                 iterations,
                                 integration_tol,
                                 integration_tol);
  }

  // Used for integrating assimilation over a range of depths within the plant
  quadrature::QAG integrator;

};

} // namespace plant

#endif
