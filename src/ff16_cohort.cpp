#include <plant.h>

// Helpers for FF16 model

// [[Rcpp::export]]
double cohort_schedule_max_time_default__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment>& p) {
   return plant::cohort_schedule_max_time_default<plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment>& p) {
   return plant::cohort_schedule_default<plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment>& p) {
   return plant::make_cohort_schedule<plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment> >(p);
}
