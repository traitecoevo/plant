// -*-c++-*-
#ifndef PLANT_PLANT_SCM_UTILS_H_
#define PLANT_PLANT_SCM_UTILS_H_

#include <plant/cohort_schedule.h>
#include <plant/disturbance.h>

namespace plant {

std::vector<double> cohort_schedule_times_default(double max_time);

template <typename Parameters>
double cohort_schedule_max_time_default(const Parameters& p) {
  Disturbance d(p.disturbance_mean_interval);
  return d.cdf(p.control.schedule_patch_survival);
}

template <typename Parameters>
CohortSchedule cohort_schedule_default(const Parameters& p) {
  const double max_time = cohort_schedule_max_time_default(p);
  CohortSchedule schedule(0);
  schedule.r_set_max_time(max_time);
  std::vector<double> times = cohort_schedule_times_default(max_time);
  if (times.size() < 1) {
    util::stop("Did not generate any times, surprisingly");
  }
  return schedule.expand(p.size(), times);
}

template <typename Parameters>
CohortSchedule make_cohort_schedule(const Parameters& p) {
  CohortSchedule ret(p.size());
  ret.r_set_max_time(p.cohort_schedule_max_time);
  ret.set_times(p.cohort_schedule_times);
  ret.r_set_ode_times(p.cohort_schedule_ode_times);
  return ret;
}

}

#endif
