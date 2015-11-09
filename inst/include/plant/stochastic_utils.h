// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_UTILS_H_
#define PLANT_PLANT_STOCHASTIC_UTILS_H_

#include <plant/cohort_schedule.h>
#include <plant/disturbance.h>
#include <plant/scm_utils.h>

namespace plant {

template <typename Parameters>
CohortSchedule make_empty_stochastic_schedule(const Parameters& p) {
  CohortSchedule ret(p.size());
  ret.r_set_max_time(p.cohort_schedule_max_time);
  return ret;
}

}

#endif
