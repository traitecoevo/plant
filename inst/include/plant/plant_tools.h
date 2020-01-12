// -*-c++-*-
#ifndef PLANT_PLANT_TOOLS_H_
#define PLANT_PLANT_TOOLS_H_

#include <plant.h>
#include <plant/uniroot.h>

namespace plant {
namespace tools {
Environment fixed_environment(double canopy_openness,
                              double height_max=150.0);
template <typename T, typename E>
double lcp_whole_plant(Plant<T,E> p) {
  auto target = [&] (double x) mutable -> double {
    Environment env = fixed_environment(x);
    p.compute_rates(env);
    return p.net_mass_production_dt(env);
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

// // These are only here because I really want somewhere after the Rcpp
// // inclusion.
// template <typename T, typename E>
// PlantPlus<T,E> plant_to_plant_plus(Plant<T,E> p, SEXP environment) {
//   if (environment == R_NilValue) {
//     return PlantPlus<T,E>(p);
//   } else {
//     return PlantPlus<T,E>(p, Rcpp::as<plant::Environment>(environment));
//   }
// }

// template <typename T, typename E>
// Plant<T,E> plant_plus_to_plant(PlantPlus<T,E> p) {
//   return p.to_plant();
// }

}

#endif
