// -*-c++-*-
#ifndef PLANT_PLANT_SCM_UTILS_H_
#define PLANT_PLANT_SCM_UTILS_H_

#include <plant/control.h>
#include <plant/cohort_schedule.h>

namespace plant {

template <typename SpeciesParameters>
CohortSchedule cohort_schedule_default(const SpeciesParameters& p, const Control& c) {
  const double max_time = c.patch_max_lifetime;
  CohortSchedule schedule(0);
  schedule.r_set_max_time(max_time);
  std::vector<double> times = p.cohort_schedule_times_default;
  if (times.size() < 1) {
    util::stop("Did not generate any times, surprisingly. Is `species_parameter$cohort_schedule_times_default` empty?");
  }
  return schedule.expand(p.size(), times);
}

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
