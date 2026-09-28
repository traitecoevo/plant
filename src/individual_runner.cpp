
#include <plant.h>

// [[Rcpp::export]]
plant::Internals<double> FF16_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::FF16_Strategy,plant::FF16_Environment>>& obj) {
  return obj.get_system().individual.r_internals();
}

// [[Rcpp::export]]
plant::Internals<double> K93_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::K93_Strategy, plant::K93_Environment>>& obj) {
  return obj.get_system().individual.r_internals();
}

// [[Rcpp::export]]
plant::Internals<double> TF24_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::TF24_Strategy<double>, plant::TF24_Environment<double>>>& obj) {
  return obj.get_system().individual.r_internals();
}



// [[Rcpp::export]]
plant::Internals<double> TF24f_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::TF24f_Strategy<double>, plant::TF24_Environment<double>>>& obj) {
  return obj.get_system().individual.r_internals();
}
