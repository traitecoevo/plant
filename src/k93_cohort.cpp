// Built from  src/ff16r_cohort.cpp on Wed Aug 12 15:33:08 2020 using the scaffolder, from the strategy:  FF16r
#include <plant.h>

// Helpers for K93 model

// [[Rcpp::export]]
double cohort_schedule_max_time_default__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Strategy,plant::K93_Environment>& p) {
   return plant::cohort_schedule_max_time_default<plant::Parameters<plant::K93_Strategy,plant::K93_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Strategy,plant::K93_Environment>& p) {
   return plant::cohort_schedule_default<plant::Parameters<plant::K93_Strategy,plant::K93_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___K93__K93_Env(const plant::Parameters<plant::K93_Strategy,plant::K93_Environment>& p) {
   return plant::make_cohort_schedule<plant::Parameters<plant::K93_Strategy,plant::K93_Environment> >(p);
}
