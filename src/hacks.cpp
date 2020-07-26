#include <plant.h>

// From RcppR6 

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

// FF16r Strategy ---------------------------------------------

// [[Rcpp::export]]
double cohort_schedule_max_time_default__Parameters___FF16r__FF16r_Env(const plant::Parameters<plant::FF16r_Strategy,plant::FF16r_Environment>& p) {
   return plant::cohort_schedule_max_time_default<plant::Parameters<plant::FF16r_Strategy,plant::FF16r_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___FF16r__FF16r_Env(const plant::Parameters<plant::FF16r_Strategy,plant::FF16r_Environment>& p) {
   return plant::cohort_schedule_default<plant::Parameters<plant::FF16r_Strategy,plant::FF16r_Environment> >(p);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___FF16r__FF16r_Env(const plant::Parameters<plant::FF16r_Strategy,plant::FF16r_Environment>& p) {
   return plant::make_cohort_schedule<plant::Parameters<plant::FF16r_Strategy,plant::FF16r_Environment> >(p);
}

// K93 Strategy ---------------------------------------------
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
