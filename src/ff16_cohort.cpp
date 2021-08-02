#include <plant.h>

// Helpers for FF16 model

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___FF16__FF16_Env(
    const plant::SpeciesParameters<plant::FF16_Strategy> &p,
    const plant::Control &c) {
  return plant::cohort_schedule_default<
      plant::SpeciesParameters<plant::FF16_Strategy>>(p, c);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___FF16__FF16_Env(
    const plant::SpeciesParameters<plant::FF16_Strategy> &p,
    const plant::Control &c) {
  return plant::make_cohort_schedule<
      plant::SpeciesParameters<plant::FF16_Strategy>>(p, c);
}
