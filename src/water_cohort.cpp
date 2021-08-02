// Built from  src/ff16_cohort.cpp on Mon Jul 19 11:01:04 2021 using the
// scaffolder, from the strategy:  FF16
#include <plant.h>

// Helpers for Water model

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___Water__FF16_Env(
    const plant::SpeciesParameters<plant::Water_Strategy> &p,
    const plant::Control &c) {
  return plant::cohort_schedule_default<
      plant::SpeciesParameters<plant::Water_Strategy>>(p, c);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___Water__FF16_Env(
    const plant::SpeciesParameters<plant::Water_Strategy> &p,
    const plant::Control &c) {
    return plant::make_cohort_schedule<
        plant::SpeciesParameters<plant::Water_Strategy>>(p, c);
  }
