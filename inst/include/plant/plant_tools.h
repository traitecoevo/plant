// -*-c++-*-
#ifndef PLANT_PLANT_TOOLS_H_
#define PLANT_PLANT_TOOLS_H_

#include <plant.h>

namespace plant {
namespace tools {
Environment fixed_environment(double canopy_openness,
                              double height_max=150.0);
double lcp_whole_plant(PlantPlus<FFW16_Strategy> p);
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
