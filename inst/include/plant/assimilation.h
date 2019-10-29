// -*-c++-*-
#ifndef PLANT_PLANT_STRATEGY_H_
#define PLANT_PLANT_STRATEGY_H_

#include <memory>
#include <plant/control.h>
#include <plant/qag_internals.h> // quadrature::intervals_type
#include <plant/internals.h> // quadrature::intervals_type

namespace plant {


template <typename T> class Assimilation {
  public:
  typedef T strategy_type;

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  double assimilation(const Environment& environment, double height,
                      double area_leaf, bool reuse_intervals);

  double compute_assimilation(double z, double height, const strategy_type strategy, 
                              const Environment& environment) const;
  
  // [eqn 12] Gross annual CO2 assimilation
  //
  // NOTE: In contrast with Daniel's implementation (but following
  // Falster 2012), we do not normalise by a_y*a_bio here.
  double assimilation(const Environment& environment,
                      double height,
                      double area_leaf,
                      bool reuse_intervals) {
    const bool over_distribution = control.plant_assimilation_over_distribution;
    const double x_min = 0.0, x_max = over_distribution ? 1.0 : height;

    double A = 0.0;

    if (control.plant_assimilation_adaptive && reuse_intervals) {
      A = control.integrator.integrate_with_last_intervals(f, x_min, x_max);
    } else {
      A = control.integrator.integrate(compute_assimilation, x_min, x_max);
    }

    return area_leaf * A;
  }

}

template <typename T> class AssimilationOverDistribution: public Assimilation {
  public:
  typedef T strategy_type;

  double compute_assimilation(strategy_type strategy, double z, double height,
                              const Environment& environment) const {
    return strategy.assimilation_leaf(environment.canopy_openness(z)) * q(z, height);
  }
}

template <typename T> class AssimilationSingleLayer: public Assimilation {
  public:
  typedef T strategy_type;

  double compute_assimilation(strategy_type strategy, double p, double height,
                              const Environment& environment) const {
    return strategy.assimilation_leaf(environment.canopy_openness(Qp(p, height)));
  }
}

}
