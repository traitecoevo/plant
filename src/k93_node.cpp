#include <plant.h>

// Helpers for K93 model

// [[Rcpp::export]]
plant::Schedule schedule_default__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Model,plant::K93_Environment>& p) {
   return plant::schedule_default<plant::Parameters<plant::K93_Model,plant::K93_Environment> >(p);
}

// [[Rcpp::export]]
plant::Schedule make_schedule__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Model,plant::K93_Environment>& p) {
   return plant::make_schedule<plant::Parameters<plant::K93_Model,plant::K93_Environment> >(p);
}
