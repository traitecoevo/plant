// -*-c++-*-
#ifndef PLANT_PLANT_ES20_ASSIMILATION_H_
#define PLANT_PLANT_ES20_ASSIMILATION_H_

#include <memory>
#include <plant/models/assimilation.h>
#include <plant/models/es20_environment.h>
#include <plant/control.h>
#include <plant/qag.h> // quadrature::intervals_type

#include <iostream>

namespace plant {

class ES20_Assimilation: public Assimilation<ES20_Environment> {
public:

  // ES20_Assimilation() : a_p1(0.0), a_p2(0.0), eta(0.0) { };
  // ES20_Assimilation(double a_p1, double a_p2, double eta) : a_p1(a_p1), a_p2(a_p2), eta(eta) { };

  // [eqn 12] Gross annual CO2 assimilation

  // NOTE: In contrast with Daniel's implementation (but following
  // Falster 2012), we do not normalise by a_y*a_bio here.
  double assimilate(Control& control,
                    const ES20_Environment& environment,
                    double height,
                    double area_leaf,
                    bool reuse_intervals) {
    const bool over_distribution = control.plant_assimilation_over_distribution;
    const double x_min = 0.0, x_max = over_distribution ? 1.0 : height;

    double A = 0.0;

    // Find the stress
    double stress = environment.getStress();
    
    // std::cout << "TIME:   " << environment.time << std::endl;
    // std::cout << "leaf area:  " << area_leaf << "  height:   " << height << "  stress:    " << stress << std::endl;



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

    return area_leaf * A * stress;
  }

};

} // namespace plant

#endif
