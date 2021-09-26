// Built from  src/ff16_cohort.cpp on Mon Jul 19 11:01:04 2021 using the scaffolder, from the strategy:  FF16
#include <plant.h>

// Helpers for FF16w model

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___FF16w__FF16_Env(const plant::Parameters<plant::FF16w_Strategy,plant::FF16_Environment>& p) {
   return plant::cohort_schedule_default<plant::Parameters<plant::FF16w_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___FF16w__FF16_Env(const plant::Parameters<plant::FF16w_Strategy,plant::FF16_Environment>& p) {
   return plant::make_cohort_schedule<plant::Parameters<plant::FF16w_Strategy,plant::FF16_Environment> >(p);
}
