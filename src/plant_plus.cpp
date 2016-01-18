#include <plant/plant_plus.h>
#include <plant/ff16_strategy.h>
#include <plant/ff16r_strategy.h>
#include <plant/plant_tools.h>

// TODO: These will disappear once the templated functions are
// implemented (see RcppR6 #23, but also #22 which would be required
// to make this nice to use).
//
// TODO: This actually can't easily be done with the RcppR6 approach
// because the different strategies aren't a class that is templated;
// instead we'll need to manually list classes (added in #23).

// [[Rcpp::export]]
plant::PlantPlus<plant::FF16_Strategy>
FF16_plant_to_plant_plus(plant::Plant<plant::FF16_Strategy> p,
                          SEXP environment) {
  return plant::plant_to_plant_plus(p, environment);
}

// [[Rcpp::export]]
plant::PlantPlus<plant::FF16r_Strategy>
FF16r_plant_to_plant_plus(plant::Plant<plant::FF16r_Strategy> p,
                          SEXP environment) {
  return plant::plant_to_plant_plus(p, environment);
}
