// Built from  src/ff16_cohort.cpp on Thu Oct 29 11:14:41 2020 using the scaffolder, from the strategy:  FF16
#include <plant.h>

// Helpers for FF16ppa model

// [[Rcpp::export]]
double cohort_schedule_max_time_default__Parameters___FF16ppa__FF16_Env(const plant::Parameters<plant::FF16ppa_Strategy,plant::FF16_Environment>& p) {
   return plant::cohort_schedule_max_time_default<plant::Parameters<plant::FF16ppa_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___FF16ppa__FF16_Env(const plant::Parameters<plant::FF16ppa_Strategy,plant::FF16_Environment>& p) {
   return plant::cohort_schedule_default<plant::Parameters<plant::FF16ppa_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___FF16ppa__FF16_Env(const plant::Parameters<plant::FF16ppa_Strategy,plant::FF16_Environment>& p) {
   return plant::make_cohort_schedule<plant::Parameters<plant::FF16ppa_Strategy,plant::FF16_Environment> >(p);
}
