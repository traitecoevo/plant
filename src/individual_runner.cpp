
#include <plant.h>

// [[Rcpp::export]]
plant::Internals<double> FF16_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::FF16_Model<double>,plant::FF16_Environment<double>>>& obj) {
  return obj.get_system().individual.r_internals();
}

// [[Rcpp::export]]
plant::Internals<double> K93_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::K93_Model, plant::K93_Environment>>& obj) {
  return obj.get_system().individual.r_internals();
}

// [[Rcpp::export]]
plant::Internals<double> TF24_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::TF24_Model<double>, plant::TF24_Environment<double>>>& obj) {
  return obj.get_system().individual.r_internals();
}



// [[Rcpp::export]]
plant::Internals<double> TF24f_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::TF24f_Model<double>, plant::TF24_Environment<double>>>& obj) {
  return obj.get_system().individual.r_internals();
}

// [[Rcpp::export]]
plant::Internals<double> TF24t_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::TF24t_Model<double>, plant::TF24_Environment<double>>>& obj) {
  return obj.get_system().individual.r_internals();
}

// [[Rcpp::export]]
plant::Internals<double> TF24w_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::TF24w_Model<double>, plant::TF24_Environment<double>>>& obj) {
  return obj.get_system().individual.r_internals();
}
