// -*-c++-*-
#ifndef PLANT_PLANT_SCM_UTILS_H_
#define PLANT_PLANT_SCM_UTILS_H_

#include <plant/node_schedule.h>

namespace plant {

std::vector<double> node_schedule_times_default(double max_time);

NodeSchedule make_mutant_schedule(const size_t n_mutants, const double max_time, 
                                  const std::vector<double> step_history);

template <typename Parameters>
NodeSchedule node_schedule_default(const Parameters& p) {
  const double max_time = p.max_patch_lifetime;
  NodeSchedule schedule(0);
  schedule.r_set_max_time(max_time);
  std::vector<double> times = node_schedule_times_default(max_time);
  if (times.size() < 1) {
    util::stop("Did not generate any times, surprisingly");
  }
  return schedule.expand(p.size(), times);
}

template <typename Parameters>
NodeSchedule make_node_schedule(const Parameters& p) {
  NodeSchedule ret(p.size());
  ret.r_set_max_time(p.max_patch_lifetime);
  ret.set_times(p.node_schedule_times);
  ret.r_set_ode_times(p.node_schedule_ode_times);
  return ret;
}

}

#endif
