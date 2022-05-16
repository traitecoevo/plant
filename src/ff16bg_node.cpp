// Built from  src/ff16_cohort.cpp on Fri Oct 30 11:43:30 2020 using the scaffolder, from the strategy:  FF16
#include <plant.h>

// Helpers for FF16bg model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___FF16bg__FF16_Env(const plant::Parameters<plant::FF16bg_Strategy,plant::FF16_Environment>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::FF16bg_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___FF16bg__FF16_Env(const plant::Parameters<plant::FF16bg_Strategy,plant::FF16_Environment>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::FF16bg_Strategy,plant::FF16_Environment> >(p);
}
