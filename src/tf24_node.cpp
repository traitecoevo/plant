#include <plant.h>

// Helpers for TF24 model

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___TF24__TF24_Env(const plant::Parameters<plant::TF24_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::schedule_default<plant::Parameters<plant::TF24_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___TF24__TF24_Env(const plant::Parameters<plant::TF24_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::make_schedule<plant::Parameters<plant::TF24_Model<double>,plant::TF24_Environment<double>> >(p);
}
