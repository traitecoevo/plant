// -*-c++-*-
#ifndef PLANT_PLANT_SCM_UTILS_H_
#define PLANT_PLANT_SCM_UTILS_H_

#include <plant/control.h>
#include <plant/cohort_schedule.h>

namespace plant {

template <typename SpeciesParameters>
CohortSchedule make_cohort_schedule(const SpeciesParameters& p, const Control& c) {
  CohortSchedule ret(p.size());
  ret.r_set_max_time(c.patch_max_lifetime);
  ret.set_times(p.cohort_schedule_times);
  ret.r_set_ode_times(p.cohort_schedule_ode_times);
  return ret;
}

}

#endif
