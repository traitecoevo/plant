// Built from src/tf24_node.cpp on Thu Jun 25 06:01:09 2026 using the scaffolder, from the strategy: TF24
// Built from  src/ff16_node.cpp on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
#include <plant.h>

// Helpers for TF24f model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___TF24f__TF24_Env(const plant::Parameters<plant::TF24f_Strategy,plant::TF24_Environment>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::TF24f_Strategy,plant::TF24_Environment> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___TF24f__TF24_Env(const plant::Parameters<plant::TF24f_Strategy,plant::TF24_Environment>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::TF24f_Strategy,plant::TF24_Environment> >(p);
}
