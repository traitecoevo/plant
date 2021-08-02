// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_UTILS_H_
#define PLANT_PLANT_STOCHASTIC_UTILS_H_

#include <plant/cohort_schedule.h>
#include <plant/scm_utils.h>

namespace plant {

template <typename SpeciesParameters>
CohortSchedule make_empty_stochastic_schedule(const SpeciesParameters& p, const Control& c) {
  CohortSchedule ret(p.size());
  ret.r_set_max_time(c.patch_max_lifetime);
  return ret;
}

}

#endif
