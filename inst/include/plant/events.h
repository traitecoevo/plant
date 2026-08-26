// -*-c++-*-
#ifndef PLANT_PLANT_EVENTS_H_
#define PLANT_PLANT_EVENTS_H_

#include <plant/node_schedule.h>
#include <plant/util.h>

#include <string>
#include <vector>

namespace plant {

// Name <-> tag for the R-facing event types (issue #522). The strings are the
// stable interface: they appear in `events()` output and in saved runs, so
// treat them as an API and add to them rather than renaming them.
EventType event_type_from_string(const std::string& name);
std::string event_type_to_string(EventType type);

// How many parameters a type expects, so a malformed event is rejected where
// the user can still see which one it was rather than deep inside a run.
size_t event_type_n_params(EventType type);

// The R-facing description of a run's discrete events.
//
// Parallel vectors rather than a vector of structs, because this crosses the R
// boundary as an RcppR6 `list:` class: plain data exposed once, rather than one
// binding per (strategy, environment) pair, which is what a queue of
// polymorphic actions over the templated Patch would have cost.
//
// This is a wire format only -- the schedule is the single source of truth
// during a run. Events goes in when the SCM is built and comes back out of the
// schedule afterwards, so a schedule refined mid-run round-trips.
class Events {
public:
  Events() {}

  size_t size() const { return time.size(); }
  // Checked on every crossing from R (RcppR6 validator), so keep it cheap.
  void validate();

  std::vector<double> time;
  std::vector<std::string> type;
  // 1-based, matching R. Meaningful only for node introductions; every other
  // type acts on the patch as a whole and ignores it.
  std::vector<size_t> species_index;
  // Per-type payload; see event_type_n_params().
  std::vector<std::vector<double> > params;
};

// Events -> queue entries. Validates species indices against the number of
// species, which Events itself cannot know.
std::vector<NodeScheduleEvent> to_schedule_events(const Events& events,
                                                  size_t n_species);

// Queue entries -> Events, for reading a (possibly refined) schedule back out.
Events events_from_schedule_events(const std::vector<NodeScheduleEvent>& events);

}

#endif
