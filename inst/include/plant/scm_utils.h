// -*-c++-*-
#ifndef PLANT_PLANT_SCM_UTILS_H_
#define PLANT_PLANT_SCM_UTILS_H_

#include <plant/node_schedule.h>
#include <plant/events.h>

namespace plant {

std::vector<double> node_schedule_times_default(double max_time);

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
  ret.r_set_ode_times(p.ode_times);
  return ret;
}

// As above, but with the schedule taken from an explicit event list rather
// than from p.node_schedule_times. An empty list means "no events supplied",
// and falls back to the parameters -- which is what keeps every existing
// caller on exactly the path it was on before events existed (#522).
template <typename Parameters>
NodeSchedule make_node_schedule(const Parameters& p, const Events& events) {
  if (events.size() == 0) {
    return make_node_schedule(p);
  }
  validate_event_horizon(events, p.max_patch_lifetime);
  NodeSchedule ret(p.size());
  ret.r_set_max_time(p.max_patch_lifetime);
  ret.set_all_events(to_schedule_events(events, p.size()));
  ret.r_set_ode_times(p.ode_times);
  return ret;
}

}

#endif
