// -*-c++-*-
#ifndef PLANT_PLANT_NODE_SCHEDULE_H_
#define PLANT_PLANT_NODE_SCHEDULE_H_

#include <RcppCommon.h> // SEXP
#include <plant/util.h>

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
};

class NodeSchedule {
public:
  typedef NodeScheduleEvent Event;
  NodeSchedule(size_t n_species_);
  size_t size() const;
  size_t get_n_species() const;
  NodeSchedule expand(size_t n_extra, std::vector<double> times);
  void clear_times(size_t species_index);
  void set_times(const std::vector<double>& times_, size_t species_index);
  void set_times(const std::vector<std::vector<double> >& times);
  std::vector<double> times(size_t species_index) const;
  // Whole-queue access, in schedule order, for round-tripping through the
  // R-facing `Events` wire format (see events.h). set_all_events() replaces
  // every event, introductions included, so it is the one entry point that
  // does not go through the per-species times interface.
  std::vector<Event> get_events() const;
  void set_all_events(const std::vector<Event>& events_);
  void reset();
  void pop();
  Event next_event() const;
  size_t remaining() const;

  double get_max_time() const;
  std::vector<std::vector<double> > get_times() const;
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
  NodeSchedule r_copy() const;

private:
  typedef std::list<Event>::iterator events_iterator;
  typedef std::list<Event>::const_iterator events_const_iterator;

  events_iterator add_time(double times, size_t species_index,
                           events_iterator it);
  events_iterator insert_event(const Event& e);
  void distribute_ode_times();

  size_t n_species;
  std::list<Event> events;
  std::list<Event> queue;
  double max_time;
  std::vector<double> ode_times;
  bool use_ode_times;
};

}

#endif
