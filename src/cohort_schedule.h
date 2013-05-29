// -*-c++-*-
#ifndef TREE_COHORT_SCHEDULE_H_
#define TREE_COHORT_SCHEDULE_H_

#include <Rcpp.h>

// NOTE: This needs copy/assignment constructors just because we need
// to copy the "next" pointer (though at present this i just done.  An
// alternative way of doing this would be for reset() to dump things
// into a Queue that we pop off?
namespace model {

class CohortSchedule {
public:
  class Event;
  CohortSchedule(size_t n_cohort_types);
  size_t size() const;
  size_t types() const;

  void clear_times(size_t cohort_index);
  void set_times(std::vector<double> times, size_t cohort_index);
  std::vector<double> times(size_t cohort_index) const;
  void reset();
  void pop();
  Event next_event() const;
  double next_time() const;
  size_t remaining() const;

  // * R interface:
  void r_clear_times(size_t cohort_index);
  void r_set_times(std::vector<double> times, size_t cohort_index);
  std::vector<double> r_times(size_t cohort_index) const;

private:
  typedef std::list<Event>::iterator events_iterator;
  typedef std::list<Event>::const_iterator events_const_iterator;

  events_iterator add_time(double times, size_t cohort_index,
			   events_iterator it);

  size_t n_cohort_types;
  std::list<Event> events;
  std::list<Event> queue;
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
