
#include <plant.h>

// Technical debt: (See RcppR6 #23 and plant #164)
// [[Rcpp::export]]
plant::Internals FF16_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::PlantRunner<plant::FF16_Strategy,plant::FF16_Environment>>& obj) {
  return obj.obj.plant.r_internals();
}

// [[Rcpp::export]]
plant::Internals FF16r_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::PlantRunner<plant::FF16r_Strategy, plant::FF16r_Environment>>& obj) {
  return obj.obj.plant.r_internals();
}


