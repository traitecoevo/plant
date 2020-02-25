
#include <plant.h>

// Technical debt: (See RcppR6 #23 and plant #164)
// [[Rcpp::export]]
plant::Internals FF16_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::PlantRunner<plant::FF16_Strategy,plant::LightEnvironment>>& obj) {
  return obj.obj.plant.r_internals();
}
