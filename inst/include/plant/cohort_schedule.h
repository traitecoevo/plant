// -*-c++-*-
#ifndef PLANT_PLANT_COHORT_SCHEDULE_H_
#define PLANT_PLANT_COHORT_SCHEDULE_H_

#include <RcppCommon.h> // SEXP
#include <plant/util.h>

// The "times" methods (set_times, times) refer to the *introduction*
// times.  As such, this really needs that at least one species has an
// introduction time of zero.  This is only enforced in that the model
// will just refuse to run if started with the incorrect time.

namespace plant {

// This could be done via a list object, but I think this is OK for
// now.  The main reason for keeping this as a separate class is it
// only makes sense to have a nontrivial constructor, and that's not
// yet supported for RcppR6 lists.
class CohortScheduleEvent {
public:
  CohortScheduleEvent(double introduction, size_t species_index_)
      : species_index(species_index_) {
    times.push_back(introduction);
  }
  size_t species_index_raw() const { return species_index; }
  double time_introduction() const { return times.front(); }
  double time_end() const { return times.back(); }

  size_t species_index;
  std::vector<double> times;
};

class CohortSchedule {
public:
  typedef CohortScheduleEvent Event;
  CohortSchedule(size_t n_species_);
  size_t size() const;
  size_t get_n_species() const;

  CohortSchedule expand(size_t n_extra, std::vector<double> times);

  void clear_times(size_t species_index);
  void set_times(const std::vector<double> &times_, size_t species_index);
  void set_times(const std::vector<std::vector<double>> &times);
  std::vector<double> times(size_t species_index) const;
  void reset();
  void pop();
  Event next_event() const;
  size_t remaining() const;

  double get_max_time() const;
  std::vector<std::vector<double>> get_times() const;
  bool using_ode_times() const;

  // * R interface:
  void r_clear_times(util::index species_index);
  std::vector<double> r_times(util::index species_index) const;
  void r_set_times(std::vector<double> times_, util::index species_index);
  void r_set_max_time(double x);
  std::vector<double> r_ode_times() const;
  void r_set_ode_times(std::vector<double> x);
  void r_clear_ode_times();
  void r_set_use_ode_times(bool x);
  SEXP r_all_times() const;
  void r_set_all_times(SEXP x);
  CohortSchedule r_copy() const;

private:
  typedef std::list<Event>::iterator events_iterator;
  typedef std::list<Event>::const_iterator events_const_iterator;

  events_iterator add_time(double times, size_t species_index,
                           events_iterator it);
  void distribute_ode_times();

  size_t n_species;
  std::list<Event> events;
  std::list<Event> queue;
  double max_time;
  std::vector<double> ode_times;
  bool use_ode_times;
};

} // namespace plant

#endif
