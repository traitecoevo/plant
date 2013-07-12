#include "cohort_schedule.h"

#include "util.h"

namespace model {

CohortSchedule::CohortSchedule(size_t n_species_)
  : n_species(n_species_),
    max_time(R_PosInf) {
  reset();
}

size_t CohortSchedule::size() const {
  return events.size();
}

size_t CohortSchedule::get_n_species() const {
  return n_species;
}

// NOTE: See note in cohort_schedule.h:CohortSchedule::Event for why
// the cast is needed here.
void CohortSchedule::clear_times(size_t species_index) {
  events_iterator e = events.begin();
  while (e != events.end()) {
    if (e->species_index == species_index)
      e = events.erase(e);
    else
      ++e;
  }
  reset();
}

void CohortSchedule::set_times(std::vector<double> times_,
			       size_t species_index) {
  clear_times(species_index);
  events_iterator e = events.begin();
  for (std::vector<double>::const_iterator t = times_.begin();
       t != times_.end(); ++t)
    e = add_time(*t, species_index, e);
  reset();
}

std::vector<double> CohortSchedule::times(size_t species_index) const {
  std::vector<double> ret;
  for (events_const_iterator e = events.begin(); e != events.end(); ++e)
    if (e->species_index == species_index)
      ret.push_back(e->time_introduction());
  return ret;
}

void CohortSchedule::reset() {
  queue = events;
  events_iterator e = queue.begin();
  while (e != queue.end()) {
    events_iterator e_next = e;
    e_next++;
    e->times.push_back(e_next == queue.end() ?
		       max_time : e_next->time_introduction());
    e = e_next;
  }
}

void CohortSchedule::pop() {
  if (queue.empty())
    ::Rf_error("Attempt to pop empty queue");
  queue.pop_front();
}

CohortSchedule::Event CohortSchedule::next_event() const {
  if (queue.empty())
    ::Rf_error("All events completed");
  return queue.front();
}

size_t CohortSchedule::remaining() const {
  return queue.size();
}

// * R interface
void CohortSchedule::r_clear_times(size_t species_index) {
  clear_times(util::check_bounds_r(species_index, n_species));
}

void CohortSchedule::r_set_times(std::vector<double> times_,
				 size_t species_index) {
  if (!util::is_sorted(times_.begin(), times_.end()))
    ::Rf_error("Times must be sorted (increasing)");
  if (times_.front() < 0)
    ::Rf_error("First time must nonnegative");
  if (times_.back() > max_time)
    ::Rf_error("Times cannot be greater than max_time");
  set_times(times_, util::check_bounds_r(species_index, n_species));
}

std::vector<double> CohortSchedule::r_times(size_t species_index) const {
  return times(util::check_bounds_r(species_index, n_species));
}

double CohortSchedule::r_max_time() const {
  return max_time;
}

void CohortSchedule::r_set_max_time(double x) {
  if (x < 0)
    ::Rf_error("max_time must be nonnegative");
  if (x < events.back().time_introduction())
    ::Rf_error("max_time must be at least the final scheduled time");
  max_time = x;
}

// * Private methods
CohortSchedule::events_iterator
CohortSchedule::add_time(double time, size_t species_index,
			 events_iterator it) {
  Event e(time, species_index);
  it = events.begin();
  while (it != events.end() && time > it->time_introduction())
    it++;
  it = events.insert(it, e);
  return it;
}

}
