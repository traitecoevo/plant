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

// NOTE: The r_cohort()/r_set_cohort() functions are here because I
// want to use NA_INTEGER for an impossible cohort.  Rcpp converts
// size_t -> numeric, though NA_REAL won't stick as a NA value there
// either.  This approach isolates the ugly casts to this class, and
// means if I change the behaviour I can just do it here.  Note we
// also add 1 so that we move to base-1 indexing.
class CohortSchedule::Event {
public:
  Event(double time_, size_t cohort_) : time(time_), cohort(cohort_) {}
  static Event blank(double time_) {
    Event e(time_, static_cast<size_t>(NA_INTEGER));
    return e;
  }
  int r_cohort() const {
    int ret = static_cast<int>(cohort);
    if (ret != NA_INTEGER)
      ret++;
    return ret;
  }
  void r_set_cohort(int x) {cohort = static_cast<size_t>(x);}

  // Better than providing an index to a cohort could be to provide an
  // iterator to the underlying cohorts.  This could be templated to
  // provide a list of different times for different things?
  double time;
  size_t cohort;
};

}

RCPP_EXPOSED_CLASS(model::CohortSchedule)
RCPP_EXPOSED_CLASS(model::CohortSchedule::Event)

#endif
