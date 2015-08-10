// -*-c++-*-
#ifndef PLANT_PLANT_TOOLS_H_
#define PLANT_PLANT_TOOLS_H_

#include <plant.h>
#include <plant/uniroot.h>

namespace plant {
namespace tools {
Environment fixed_environment(double canopy_openness,
                              double height_max=150.0);
template <typename T>
double lcp_whole_plant(PlantPlus<T> p) {
  auto target = [&] (double x) mutable -> double {
    Environment env = fixed_environment(x);
    p.compute_vars_phys(env);
    return p.net_mass_production_dt();
  };

  const double f1 = target(1.0);
  if (f1 < 0.0) {
    return NA_REAL;
  } else {
    const double tol = p.control().plant_seed_tol;
    const size_t max_iterations = p.control().plant_seed_iterations;
    return util::uniroot(target, 0.0, 1.0, tol, max_iterations);
  }
}
}

// These are only here because I really want somewhere after the Rcpp
// inclusion.
template <typename T>
PlantPlus<T> plant_to_plant_plus(Plant<T> p, SEXP environment) {
  if (environment == R_NilValue) {
    return PlantPlus<T>(p);
  } else {
    return PlantPlus<T>(p, Rcpp::as<plant::Environment>(environment));
  }
}

template <typename T>
Plant<T> plant_plus_to_plant(PlantPlus<T> p) {
  return p.to_plant();
}

}

#endif
