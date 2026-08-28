// -*-c++-*-
#ifndef PLANT_PLANT_EVENTS_H_
#define PLANT_PLANT_EVENTS_H_

#include <plant/node_schedule.h>
#include <plant/util.h>

#include <string>
#include <vector>

namespace plant {

// Name <-> tag for the R-facing event types (issue #628). The strings are the
// stable interface: they appear in `events()` output and in saved runs, so
// treat them as an API and add to them rather than renaming them.
EventType event_type_from_string(const std::string& name);
std::string event_type_to_string(EventType type);

EventTarget event_target_from_string(const std::string& name);
std::string event_target_to_string(EventTarget target);

// How many parameters a type expects, so a malformed event is rejected where
// the user can still see which one it was rather than deep inside a run.
size_t event_type_n_params(EventType type);

// The target a type acts on when none is named, and whether it will accept
// being narrowed to a single species. Harvest and climate extremes will (act on
// one species, or on the whole patch); a resource pulse will not, because a
// resource pool belongs to the environment rather than to any one species.
EventTarget event_type_default_target(EventType type);
bool event_type_accepts_species(EventType type);
// Whether a type will act on this target at all, and the list it will take.
bool event_type_allows_target(EventType type, EventTarget target);
std::string event_type_target_list(EventType type);


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
  // What each event acts on: "patch", "environment" or "species".
  std::vector<std::string> target;
  // 1-based, matching R. Which species, or which resource, depending on the
  // target; unread when the target is the whole patch.
  std::vector<size_t> target_index;
  // Per-type payload; see event_type_n_params().
  std::vector<std::vector<double> > params;
};

// What a run's events actually did, as opposed to what was asked of them
// (issue #628).
//
// The two are routinely different and the difference is the interesting part:
// a resource pulse is capped at what the pool can hold, so the amount that
// lands is often less than the amount requested, and harvesting a size class
// removes whatever was in that class rather than a fixed number.
// Without this the shortfall is only inferable from an accumulator, which is
// no way to answer "what did this run do".
class EventLog {
public:
  EventLog() {}
  size_t size() const { return time.size(); }
  void validate() {}

  std::vector<double> time;
  std::vector<std::string> type;
  std::vector<std::string> target;
  std::vector<size_t> target_index;
  // What the event asked for: the event's own params, verbatim.
  std::vector<std::vector<double> > requested;
  // What it achieved. Per-type, and documented by the action that fills it:
  // a pulse reports {accepted, shed}; harvest and climate extremes report
  // {fraction_applied, nodes_affected, density_removed}.
  //
  // `nodes_affected` counts numerical cohorts touched, not individuals -- it is
  // bookkeeping about the discretisation. `density_removed` is the quantity
  // actually taken out, summed over those cohorts, and is the number to read
  // when asking what the intervention did.
  std::vector<std::vector<double> > applied;
};

// One applied event, as the action reports it back to the runner.
class EventRecord {
public:
  double time = 0.0;
  EventType type = EventType::NodeIntroduction;
  EventTarget target = EventTarget::Patch;
  size_t target_index = 0;
  std::vector<double> requested;
  std::vector<double> applied;
};

EventLog event_log_from_records(const std::vector<EventRecord>& records);

// Reject events scheduled past the end of the run. Separate from validate()
// because Events crosses the R boundary without knowing the run's length.
void validate_event_horizon(const Events& events, double max_time);

// Events -> queue entries. Validates species indices against the number of
// species, which Events itself cannot know.
std::vector<NodeScheduleEvent> to_schedule_events(const Events& events,
                                                  size_t n_species);

// Queue entries -> Events, for reading a (possibly refined) schedule back out.
Events events_from_schedule_events(const std::vector<NodeScheduleEvent>& events);

}

#endif
