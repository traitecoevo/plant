#include <plant.h>

// Helpers for TF24w model

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___TF24w__TF24_Env(const plant::Parameters<plant::TF24w_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::schedule_default<plant::Parameters<plant::TF24w_Model<double>,plant::TF24_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___TF24w__TF24_Env(const plant::Parameters<plant::TF24w_Model<double>,plant::TF24_Environment<double>>& p) {
   return plant::make_schedule<plant::Parameters<plant::TF24w_Model<double>,plant::TF24_Environment<double>> >(p);
}
