// Built from  src/ff16_node.cpp on Mon Jul 19 11:01:04 2021 using the scaffolder, from the strategy:  FF16
#include <plant.h>

// Helpers for FF16w model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___FF16w__FF16_Env(const plant::Parameters<plant::FF16w_Strategy,plant::FF16_Environment>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::FF16w_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___FF16w__FF16_Env(const plant::Parameters<plant::FF16w_Strategy,plant::FF16_Environment>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::FF16w_Strategy,plant::FF16_Environment> >(p);
}
