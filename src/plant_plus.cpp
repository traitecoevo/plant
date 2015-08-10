#include <plant/plant_plus.h>
#include <plant/ffw16_strategy.h>
#include <plant/ffdev_strategy.h>
#include <plant/plant_tools.h>

// TODO: These will disappear once the templated functions are
// implemented (see RcppR6 #23, but also #22 which would be required
// to make this nice to use).
//
// TODO: it would be nice to pass environment=NULL as a default here
// but I don't see how as we definitely don't want C++'s NULL!
//
// TODO: This actually can't easily be done with the RcppR6 approach
// because the different strategies aren't a class that is templated;
// instead we'll need to manually list classes (added in #23).

// [[Rcpp::export]]
plant::PlantPlus<plant::FFW16_Strategy>
FFW16_plant_to_plant_plus(plant::Plant<plant::FFW16_Strategy> p,
                          SEXP environment) {
  return plant::plant_to_plant_plus(p, environment);
}

// [[Rcpp::export]]
plant::PlantPlus<plant::FFdev_Strategy>
FFdev_plant_to_plant_plus(plant::Plant<plant::FFdev_Strategy> p,
                          SEXP environment) {
  return plant::plant_to_plant_plus(p, environment);
}
