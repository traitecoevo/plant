// -*-c++-*-
#ifndef TREE_COHORT_SCHEDULE_H_
#define TREE_COHORT_SCHEDULE_H_

#include <Rcpp.h>

namespace model {

class CohortSchedule {
public:
  class Event;
  CohortSchedule(size_t n_cohort_types);
  size_t size() const;

  void clear_times(size_t cohort);
  void set_times(std::vector<double> times, size_t cohort);
  std::vector<double> times(size_t cohort) const;
  void reset();
  Event next_event();
  double next_time() const;

private:
  void add_time(double times, size_t cohort);

  size_t n_cohort_types;
  std::list<Event> events;
  typedef std::list<Event>::iterator events_iterator;
  typedef std::list<Event>::const_iterator events_const_iterator;
  events_iterator next;
};

// NOTE: I'm using int here over size_t in cohort because I want to
// use NA_INTEGER for an impossible cohort.  Rcpp converts size_t ->
// numeric, though NA_REAL won't stick as a NA value there either.
// This causes a few casts in code.
class CohortSchedule::Event {
public:
  Event(double time, int cohort) : time(time), cohort(cohort) {}
  static Event blank() {
    Event e(R_PosInf, NA_INTEGER);
    return e;
  }
  // Better than providing an index to a cohort could be to provide an
  // iterator to the underlying cohorts.  This could be templated to
  // provide a list of different times for different things?
  double time;
  int cohort;
};

}

RCPP_EXPOSED_CLASS(model::CohortSchedule)
RCPP_EXPOSED_CLASS(model::CohortSchedule::Event)

#endif
