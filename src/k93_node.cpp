// Built from  src/ff16r_node.cpp on Wed Aug 12 15:33:08 2020 using the scaffolder, from the strategy:  FF16r
#include <plant.h>

// Helpers for K93 model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Strategy,plant::K93_Environment>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::K93_Strategy,plant::K93_Environment> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Strategy,plant::K93_Environment>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::K93_Strategy,plant::K93_Environment> >(p);
}
