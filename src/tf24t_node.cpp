// Built from src/tf24_node.cpp on Wed Jul 22 16:42:11 2026 using the scaffolder, from the strategy: TF24
// Built from  src/ff16_node.cpp on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
#include <plant.h>

// Helpers for TF24t model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___TF24t__TF24_Env(const plant::Parameters<plant::TF24t_Strategy,plant::TF24_Environment>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::TF24t_Strategy,plant::TF24_Environment> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___TF24t__TF24_Env(const plant::Parameters<plant::TF24t_Strategy,plant::TF24_Environment>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::TF24t_Strategy,plant::TF24_Environment> >(p);
}
