#include <plant.h>

// Helpers for FF16 model

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Model<double>,plant::FF16_Environment<double>>& p) {
   return plant::schedule_default<plant::Parameters<plant::FF16_Model<double>,plant::FF16_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Model<double>,plant::FF16_Environment<double>>& p) {
   return plant::make_schedule<plant::Parameters<plant::FF16_Model<double>,plant::FF16_Environment<double>> >(p);
}
