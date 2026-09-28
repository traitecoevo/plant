// -*-c++-*-
#ifndef PLANT_PLANT_NODE_SCHEDULE_H_
#define PLANT_PLANT_NODE_SCHEDULE_H_

#include <RcppCommon.h> // SEXP
#include <plant/util.h>
#include <odelia/ode_interface.hpp>

// The "times" methods (set_times, times) refer to the *introduction*
// times.  As such, this really needs that at least one species has an
// introduction time of zero.  This is only enforced in that the model
// will just refuse to run if started with the incorrect time.

namespace plant {

// What an event *does* when the solver reaches its time (issue #628).
//
// Every event in the schedule used to be a node introduction; the tag is what
// lets resource pulses, harvests and the rest share one queue, so that
// SCM::run_next_impl's stop/apply/resume loop serves all of them.
//
// These names are deliberately taxa- and model-agnostic. This layer is shared
// by every strategy and environment, and nothing about it is specific to
// plants, to water, or to temperature: a resource pulse is water in TF24 and
// could be anything countable in a size-structured animal model, and a climate
// extreme is heat in one model and could be cold, salinity or hypoxia in
// another. Model-specific vocabulary belongs in the model, where it is
// accurate -- TF24_Environment::add_water_pulse() is this same action under
// the name that reads correctly there.
//
// The enumerator order is also the order in which events sharing a time are
// applied, so do not reorder casually. Environment events come first (they set
// the external conditions), removals next, and node introduction last so that
// a newborn's initial conditions are computed against the post-event
// environment.
enum class EventType {
  ResourcePulse = 0,
  ClimateExtreme,
  Harvest,
  NodeIntroduction
};

// What an event acts on (issue #628). Separate from the type, because the same
// action can be aimed at different scopes: harvesting a whole patch and
// harvesting one species run the same code over a different set of nodes.
//
// Deliberately no per-cohort target. A cohort has no stable address across a
// run -- nodes are appended and never removed, and refine_schedule() changes
// how many there are -- so "cohort 7" in a schedule written before the run is
// not well defined. Selecting particular cohorts is expressed as a predicate on
// their state (a size range) in the action's parameters instead, which is both
// well defined and what the size-selective harvest use case actually asks for.
enum class EventTarget {
  Patch = 0,     // the whole patch: every species, every node
  Environment,   // the external state: resource pools and drivers
  Species        // one species, named by target_index
};

// Position of a type in the within-time application order. The enumerator
// value *is* the rank; the function exists so that call sites read as intent
// rather than as an incidental cast.
inline int event_type_rank(EventType type) {
  return static_cast<int>(type);
}

// This could be done via a list object, but I think this is OK for
// now.  The main reason for keeping this as a separate class is it
// only makes sense to have a nontrivial constructor, and that's not
// yet supported for RcppR6 lists.
class NodeScheduleEvent {
public:
  NodeScheduleEvent(double introduction, size_t target_index_,
                    EventType type_ = EventType::NodeIntroduction,
                    EventTarget target_ = EventTarget::Species,
                    std::vector<double> params_ = std::vector<double>())
    : type(type_), target(target_), target_index(target_index_),
      params(params_) {
    times.push_back(introduction);
  }
  size_t species_index_raw() const {
    return target_index;
  }
  double time_introduction() const {
    return times.front();
  }
  double time_end() const {
    return times.back();
  }
  bool is_node_introduction() const {
    return type == EventType::NodeIntroduction;
  }

  EventType type;
  EventTarget target;
  // Which one, within the target: the species for a species-targeted event,
  // the resource for an environment-targeted one. Zero and unused when the
  // target is the whole patch.
  size_t target_index;
  // Per-type payload: a pulse depth, a harvested fraction, and so on. Kept as
  // a bare vector so the queue stays plain data and crosses the R boundary
  // without a templated binding per (strategy, environment) pair.
  std::vector<double> params;
  std::vector<double> times;

// One instant, and everything that happens at it.
//
// The schedule is ascending by time with ONE row per distinct time, because a
// time naming three species is one introduction: the run introduces them
// together and computes the environment once for the set. Actions at the same
// instant are held sorted by event_type_rank and apply BEFORE the
// introductions, so a newborn's initial conditions are computed against the
// post-event environment -- an ordering made a property of this structure
// rather than re-established by the run loop at every stop.
//
// `times` keeps the [t_intro, ...extra ode times..., t_end] semantics events
// carried before; ode_steps is what a replay reads, and program_within() serves it.
};

struct schedule_entry {
  double time;
  std::vector<NodeScheduleEvent> actions;
  std::vector<size_t> species;
  std::vector<double> times;
  double time_introduction() const { return times.front(); }
  double time_end() const { return times.back(); }
};

class NodeSchedule {
public:
  typedef NodeScheduleEvent Event;
  NodeSchedule(size_t n_species_);
  // Introductions, not species-times: a time naming three species is one.
  size_t size() const;
  size_t get_n_species() const;
  NodeSchedule expand(size_t n_extra, std::vector<double> times);
  void clear_times(size_t species_index);
  void set_times(const std::vector<double>& times_, size_t species_index);
  void set_times(const std::vector<std::vector<double> >& times);
  std::vector<double> times(size_t species_index) const;
  // Puts the walk back at the first introduction. The schedule itself is not
  // consumed by a run, so this restores a position and nothing else.
  // Whole-queue access, in schedule order, for round-tripping through the
  // R-facing `Events` wire format (see events.h). set_all_events() replaces
  // every event, introductions included, so it is the one entry point that
  // does not go through the per-species times interface.
  std::vector<Event> get_events() const;
  void set_all_events(const std::vector<Event>& events_);
  void reset();

  // Where the walk is, where its interval ends, and moving past it. `next()`
  // stays valid across a pop, which only moves the position.
  const schedule_entry& next() const;
  double time_end() const;
  void pop();
  size_t remaining() const;

  // The program a replay of the interval `(start, end)` takes: the state it
  // starts from, then the recorded steps strictly inside it. Empty where the
  // recording put no step in there, and empty on a schedule holding no recording.
  //
  // Boundaries are excluded because the run stops at them anyway and the
  // interval above starts from there -- so a step at one would be taken twice.
  std::vector<odelia::ode::instruction>
  program_within(double start, double end) const;

  double get_max_time() const;
  std::vector<std::vector<double> > get_times() const;
  bool using_ode_steps() const { return !ode_steps.empty(); }

  // * R interface:
  void r_clear_times(util::index species_index);
  std::vector<double> r_times(util::index species_index) const;
  void r_set_times(std::vector<double> times_, util::index species_index);
  void r_set_max_time(double x);
  std::vector<double> r_ode_times() const;
  std::vector<double> r_ode_step_sizes() const;
  // The two halves of one recording, installed together: apart, they can be
  // paired across different runs and nothing says so.
  void r_set_ode_steps(std::vector<double> times, std::vector<double> sizes);
  void r_clear_ode_steps();
  // Where the walk is, as a list: the time, the species from one, and the time
  // its interval ends. The species are one number each, counted the way R counts
  // them, where the class this replaces reported the same index twice under two
  // names.
  SEXP r_next_introduction() const;
  SEXP r_all_times() const;
  void r_set_all_times(SEXP x);
  NodeSchedule r_copy() const;

private:
  void insert(double time, size_t species_index);
  std::vector<schedule_entry>::iterator entry_at(double time);

  size_t n_species;
  // Ascending by time, one entry per distinct time, species ascending within an
  // entry -- which is the order Patch::introduced_at rebuilds off the same times,
  // so the two agree.
  std::vector<schedule_entry> schedule;
  // How far the run has got. A position rather than a consumed copy of the
  // schedule, so what the run reads and what a caller set are one object.
  size_t at;
  double max_time;
  // The recorded ODE steps this schedule replays, when it holds any. Holding
  // them is what makes a schedule a replay, so no flag says so separately.
  std::vector<odelia::ode::instruction> ode_steps;
};

}

#endif
