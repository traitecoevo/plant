#include <plant.h>

// Helpers for FF16r model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___FF16r__FF16_Env(const plant::Parameters<plant::FF16r_Strategy,plant::FF16_Environment>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::FF16r_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___FF16r__FF16_Env(const plant::Parameters<plant::FF16r_Strategy,plant::FF16_Environment>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::FF16r_Strategy,plant::FF16_Environment> >(p);
}
