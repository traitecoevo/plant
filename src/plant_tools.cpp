
#include <plant.h>

// [[Rcpp::export]]
plant::Internals FF16_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::IndividualRunner<plant::FF16_Strategy,plant::FF16_Environment>>& obj) {
  return obj.obj.individual.r_internals();
}

// [[Rcpp::export]]
plant::Internals FF16r_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::IndividualRunner<plant::FF16r_Strategy, plant::FF16_Environment>>& obj) {
  return obj.obj.individual.r_internals();
}



// [[Rcpp::export]]
plant::Internals K93_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::IndividualRunner<plant::K93_Strategy, plant::K93_Environment>>& obj) {
  return obj.obj.individual.r_internals();
}



// [[Rcpp::export]]
plant::Internals FF16w_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::IndividualRunner<plant::FF16w_Strategy, plant::FF16_Environment>>& obj) {
  return obj.obj.individual.r_internals();
}
