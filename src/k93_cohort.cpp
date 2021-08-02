// Built from  src/ff16r_cohort.cpp on Wed Aug 12 15:33:08 2020 using the
// scaffolder, from the strategy:  FF16r
#include <plant.h>

// Helpers for K93 model

// [[Rcpp::export]]
plant::CohortSchedule cohort_schedule_default__Parameters___K93__K93_Env(
    const plant::SpeciesParameters<plant::K93_Strategy> &p,
    const plant::Control &c) {
  return plant::cohort_schedule_default<
      plant::SpeciesParameters<plant::K93_Strategy>>(p, c);
}

// [[Rcpp::export]]
plant::CohortSchedule make_cohort_schedule__Parameters___K93__K93_Env(
    const plant::SpeciesParameters<plant::K93_Strategy> &p,
    const plant::Control &c) {
  return plant::make_cohort_schedule<
      plant::SpeciesParameters<plant::K93_Strategy>>(p, c);
}
