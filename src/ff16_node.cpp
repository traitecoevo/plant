#include <plant.h>

// Helpers for FF16 model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy<double>,plant::FF16_Environment<double>>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::FF16_Strategy<double>,plant::FF16_Environment<double>> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy<double>,plant::FF16_Environment<double>>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::FF16_Strategy<double>,plant::FF16_Environment<double>> >(p);
}
