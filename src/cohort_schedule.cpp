#include "cohort_schedule.h"

#include "util.h"

namespace model {

CohortSchedule::CohortSchedule(size_t n_cohort_types)
  : n_cohort_types(n_cohort_types) {
  reset();
}

size_t CohortSchedule::size() const {
  return events.size();
}

// NOTE: See note in cohort_schedule.h:CohortSchedule::Event for why
// the cast is needed here.
void CohortSchedule::clear_times(size_t cohort) {
  const size_t cohort_base_0 = util::check_bounds_r(cohort, n_cohort_types);
  events_iterator e = events.begin();
  while (e != events.end()) {
    if (e->cohort == static_cast<int>(cohort_base_0))
      e = events.erase(e);
    else
      ++e;
  }
}

void CohortSchedule::set_times(std::vector<double> times, size_t cohort) {
  if (!util::is_sorted(times.begin(), times.end()))
    ::Rf_error("Times must be sorted (increasing)");
  clear_times(cohort);
  reset();
  const size_t cohort_base_0 = util::check_bounds_r(cohort, n_cohort_types);
  for (std::vector<double>::const_iterator it = times.begin();
       it != times.end(); it++)
    add_time(*it, cohort_base_0);
  reset();
}

std::vector<double> CohortSchedule::times(size_t cohort) const {
  const size_t cohort_base_0 = util::check_bounds_r(cohort, n_cohort_types);
  std::vector<double> ret;
  for (events_const_iterator e = events.begin(); e != events.end(); e++)
    if (e->cohort == static_cast<int>(cohort_base_0))
      ret.push_back(e->time);
  return ret;
}

void CohortSchedule::reset() {
  next = events.begin();
}

CohortSchedule::Event CohortSchedule::next_event() {
  if (next == events.end())
    return CohortSchedule::Event::blank();
  else
    return *next++;
}

double CohortSchedule::next_time() const {
  return next == events.end() ?
    CohortSchedule::Event::blank().time : next->time;
}

void CohortSchedule::add_time(double time, size_t cohort) {
  Event e(time, cohort);
  while (next != events.end() && time > next->time)
    next++;
  next = events.insert(next, e);
}

}
