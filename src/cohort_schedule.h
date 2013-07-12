// -*-c++-*-
#ifndef TREE_COHORT_SCHEDULE_H_
#define TREE_COHORT_SCHEDULE_H_

#include <Rcpp.h>
#include "util.h"

// The "times" methods (set_times, times) refer to the *introduction*
// times.  As such, this really needs that at least one species has an
// introduction time of zero.  This is only enforced in that the model
// will just refuse to run if started with the incorrect time.

namespace model {

class CohortSchedule {
public:
  class Event;
  CohortSchedule(size_t n_species_);
  size_t size() const;
  size_t get_n_species() const;

  void clear_times(size_t species_index);
  void set_times(std::vector<double> times_, size_t species_index);
  std::vector<double> times(size_t species_index) const;
  void reset();
  void pop();
  Event next_event() const;
  double next_time() const;
  size_t remaining() const;

  // * R interface:
  void r_clear_times(size_t species_index);
  void r_set_times(std::vector<double> times_, size_t species_index);
  std::vector<double> r_times(size_t species_index) const;
  double r_max_time() const;
  void r_set_max_time(double x);

private:
  typedef std::list<Event>::iterator events_iterator;
  typedef std::list<Event>::const_iterator events_const_iterator;

  events_iterator add_time(double times, size_t species_index,
			   events_iterator it);

  size_t n_species;
  std::list<Event> events;
  std::list<Event> queue;
  double max_time;
};

class CohortSchedule::Event {
public:
  Event(double introduction, size_t cohort_) : cohort(cohort_) {
    times.push_back(introduction);
  }
  static Event blank(double time_) {
    Event e(time_, util::base_1_to_0<int,size_t>(NA_INTEGER));
    return e;
  }
  int r_cohort() const {
    return util::base_0_to_1<size_t,int>(cohort);
  }
  std::vector<double> r_times() const {
    std::vector<double> ret(times.begin(), times.end());
    return ret;
  }
  double time_introduction() const {
    return times.front();
  }
  double time_end() const {
    return times.back();
  }

  size_t cohort;
  std::list<double> times;
};

}

RCPP_EXPOSED_CLASS(model::CohortSchedule)
RCPP_EXPOSED_CLASS(model::CohortSchedule::Event)

#endif
