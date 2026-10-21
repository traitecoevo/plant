#include <plant.h>

// Helpers for TF24t model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___TF24t__TF24_Env(const plant::Parameters<plant::TF24t_Strategy<double>,plant::TF24_Environment<double>>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::TF24t_Strategy<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___TF24t__TF24_Env(const plant::Parameters<plant::TF24t_Strategy<double>,plant::TF24_Environment<double>>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::TF24t_Strategy<double>,plant::TF24_Environment<double>> >(p);
}
