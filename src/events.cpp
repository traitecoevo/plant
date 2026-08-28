#include <plant/events.h>
#include <Rcpp.h>

namespace plant {

namespace {
// One table, read every way, so a name, its tag, its arity and what it may be
// aimed at cannot drift apart. Order follows EventType's, which is also the
// within-time application order.
struct TypeInfo {
  EventType type;
  const char* name;
  size_t n_params;
  EventTarget default_target;
  // Every target this type will accept. A type that acts on the abiotic state
  // cannot be aimed at a species, and one that acts on individuals cannot be
  // aimed at the environment; stating both directions here is what stops an
  // impossible combination being accepted and then quietly reinterpreted.
  bool allows_patch;
  bool allows_environment;
  bool allows_species;
};

const TypeInfo type_info[] = {
  // amount. The resource is named by target_index, and belongs to the
  // environment rather than to any one species, so a pulse cannot be narrowed
  // to one species.
  {EventType::ResourcePulse,    "resource_pulse",    1,
   EventTarget::Environment, false, true,  false},
  // intensity, duration, threshold, sensitivity
  {EventType::ClimateExtreme,   "climate_extreme",   4,
   EventTarget::Patch,       true,  false, true},
  // fraction, size_min, size_max
  {EventType::Harvest,          "harvest",           3,
   EventTarget::Patch,       true,  false, true},
  {EventType::NodeIntroduction, "node_introduction", 0,
   EventTarget::Species,     false, false, true}
};

const size_t n_type_info = sizeof(type_info) / sizeof(type_info[0]);

const TypeInfo& info_for(EventType type) {
  for (size_t i = 0; i < n_type_info; ++i) {
    if (type == type_info[i].type) {
      return type_info[i];
    }
  }
  util::stop("Unknown event type");
  return type_info[0]; // not reached
}

struct TargetName {
  EventTarget target;
  const char* name;
};

const TargetName target_names[] = {
  {EventTarget::Patch,       "patch"},
  {EventTarget::Environment, "environment"},
  {EventTarget::Species,     "species"}
};

const size_t n_target_names = sizeof(target_names) / sizeof(target_names[0]);

std::string known_type_names() {
  std::string ret;
  for (size_t i = 0; i < n_type_info; ++i) {
    ret += std::string(i > 0 ? ", " : "") + type_info[i].name;
  }
  return ret;
}
}

EventType event_type_from_string(const std::string& name) {
  for (size_t i = 0; i < n_type_info; ++i) {
    if (name == type_info[i].name) {
      return type_info[i].type;
    }
  }
  util::stop("Unknown event type '" + name + "'. Known types: " +
             known_type_names());
  return EventType::NodeIntroduction; // not reached
}

std::string event_type_to_string(EventType type) {
  return info_for(type).name;
}

size_t event_type_n_params(EventType type) {
  return info_for(type).n_params;
}

EventTarget event_type_default_target(EventType type) {
  return info_for(type).default_target;
}

bool event_type_accepts_species(EventType type) {
  return info_for(type).allows_species;
}

bool event_type_allows_target(EventType type, EventTarget target) {
  const TypeInfo& i = info_for(type);
  switch (target) {
  case EventTarget::Patch:       return i.allows_patch;
  case EventTarget::Environment: return i.allows_environment;
  case EventTarget::Species:     return i.allows_species;
  }
  return false;
}

// The targets a type will take, for an error message that says what to do
// rather than only what was wrong.
std::string event_type_target_list(EventType type) {
  const TypeInfo& i = info_for(type);
  std::string ret;
  if (i.allows_patch)       ret += std::string(ret.empty() ? "" : ", ") + "patch";
  if (i.allows_environment) ret += std::string(ret.empty() ? "" : ", ") + "environment";
  if (i.allows_species)     ret += std::string(ret.empty() ? "" : ", ") + "species";
  return ret;
}

EventTarget event_target_from_string(const std::string& name) {
  for (size_t i = 0; i < n_target_names; ++i) {
    if (name == target_names[i].name) {
      return target_names[i].target;
    }
  }
  util::stop("Unknown event target '" + name +
             "'. Known targets: patch, environment, species");
  return EventTarget::Patch; // not reached
}

std::string event_target_to_string(EventTarget target) {
  for (size_t i = 0; i < n_target_names; ++i) {
    if (target == target_names[i].target) {
      return target_names[i].name;
    }
  }
  util::stop("Unknown event target");
  return ""; // not reached
}

namespace {
// Reject a parameter the action cannot act on, here rather than mid-run. Each
// of these is silent otherwise: a NaN intensity becomes a zero-damage climate
// event, a negative sensitivity produces negative applied mortality, and an
// inverted size band harvests nothing at all -- results that look like answers.
void validate_event_params(const std::string& at, EventType type,
                           const std::vector<double>& q) {
  const auto finite = [&](double v, const char* name) {
    if (!util::is_finite(v)) {
      util::stop(at + " has a non-finite " + name);
    }
  };
  const auto nonneg = [&](double v, const char* name) {
    finite(v, name);
    if (v < 0.0) {
      util::stop(at + " has a negative " + name + " (" +
                 util::to_string(v) + ")");
    }
  };
  switch (type) {
  case EventType::ResourcePulse:
    nonneg(q.at(0), "amount");
    break;
  case EventType::ClimateExtreme:
    finite(q.at(0), "intensity");
    nonneg(q.at(1), "duration");
    finite(q.at(2), "threshold");
    nonneg(q.at(3), "sensitivity");
    break;
  case EventType::Harvest: {
    finite(q.at(0), "fraction");
    if (q.at(0) < 0.0 || q.at(0) >= 1.0) {
      util::stop(at + " has fraction " + util::to_string(q.at(0)) +
                 ", which must be in [0, 1): removing all of a cohort would "
                 "take its density to -Inf.");
    }
    nonneg(q.at(1), "size_min");
    // size_max is Inf by default, which is the whole point of the default, so
    // it is bounded rather than required finite.
    if (std::isnan(q.at(2)) || q.at(2) < q.at(1)) {
      util::stop(at + " has size_max " + util::to_string(q.at(2)) +
                 " below size_min " + util::to_string(q.at(1)) +
                 ": the band is empty, so the event would silently do nothing");
    }
    break;
  }
  case EventType::NodeIntroduction:
    break;
  }
}
}

void Events::validate() {
  const size_t n = time.size();
  if (type.size() != n || target.size() != n || target_index.size() != n ||
      params.size() != n) {
    util::stop("Events columns must all have the same length (time has " +
               util::to_string(n) + ")");
  }
  for (size_t i = 0; i < n; ++i) {
    const std::string at = "Event " + util::to_string(i + 1) + " (" +
      type[i] + ")";
    if (!util::is_finite(time[i]) || time[i] < 0.0) {
      util::stop(at + " has a non-finite or negative time");
    }
    // Also rejects an unknown name, with the list of known ones.
    const EventType t = event_type_from_string(type[i]);
    const size_t n_expected = event_type_n_params(t);
    if (params[i].size() != n_expected) {
      util::stop(at + " expects " + util::to_string(n_expected) +
                 " parameters but has " + util::to_string(params[i].size()));
    }
    const EventTarget tg = event_target_from_string(target[i]);
    if (!event_type_allows_target(t, tg)) {
      // Naming the permitted targets matters more than naming the offence: an
      // environment-aimed harvest is a modelling mistake, and the useful reply
      // is what it should have said instead.
      if (tg == EventTarget::Species) {
        util::stop(at + " cannot be aimed at a single species; it accepts " +
                   event_type_target_list(t));
      }
      if (t == EventType::NodeIntroduction) {
        util::stop(at + " must name the species being introduced");
      }
      util::stop(at + " cannot act on the " + target[i] + "; it accepts " +
                 event_type_target_list(t));
    }
    validate_event_params(at, t, params[i]);
  }
}

void validate_event_horizon(const Events& events, double max_time) {
  // An event past the end of the run is not a no-op to be dropped: it usually
  // means a time-unit slip, a truncated run, or a forcing record reused at the
  // wrong length. Silently ignoring it would give a plausible answer that is
  // missing the intervention the user asked for. Caught here, before any
  // integration, rather than surfacing later as the solver's own complaint
  // about being asked to integrate backwards.
  for (size_t i = 0; i < events.size(); ++i) {
    if (events.time[i] > max_time) {
      util::stop("Event " + util::to_string(i + 1) + " (" + events.type[i] +
                 ") occurs at time " + util::to_string(events.time[i]) +
                 ", after max_patch_lifetime = " + util::to_string(max_time) +
                 ". Events must fall within [0, max_patch_lifetime].");
    }
  }
}

std::vector<NodeScheduleEvent> to_schedule_events(const Events& events,
                                                  size_t n_species) {
  std::vector<NodeScheduleEvent> ret;
  ret.reserve(events.size());
  for (size_t i = 0; i < events.size(); ++i) {
    const EventType type = event_type_from_string(events.type[i]);
    const EventTarget target = event_target_from_string(events.target[i]);
    // 1-based on the way in, as everywhere else in the R interface.
    size_t index = 0;
    if (target != EventTarget::Patch) {
      const size_t raw = events.target_index[i];
      if (raw < 1) {
        util::stop("Event " + util::to_string(i + 1) +
                   " has target_index " + util::to_string(raw) +
                   ", which must be at least 1");
      }
      // Species can be bounds-checked here; a resource index cannot, because
      // how many an environment has is the environment's business and it is
      // not in scope. Patch::apply_event() checks it against n_resources().
      if (target == EventTarget::Species && raw > n_species) {
        util::stop("Event " + util::to_string(i + 1) +
                   " has target_index " + util::to_string(raw) +
                   ", outside 1.." + util::to_string(n_species));
      }
      index = raw - 1;
    }
    ret.push_back(NodeScheduleEvent(events.time[i], index, type, target,
                                    events.params[i]));
  }
  return ret;
}

Events events_from_schedule_events(
    const std::vector<NodeScheduleEvent>& events) {
  Events ret;
  ret.time.reserve(events.size());
  ret.type.reserve(events.size());
  ret.target.reserve(events.size());
  ret.target_index.reserve(events.size());
  ret.params.reserve(events.size());
  for (std::vector<NodeScheduleEvent>::const_iterator e = events.begin();
       e != events.end(); ++e) {
    ret.time.push_back(e->time_introduction());
    ret.type.push_back(event_type_to_string(e->type));
    ret.target.push_back(event_target_to_string(e->target));
    ret.target_index.push_back(e->target_index + 1); // back to 1-based
    ret.params.push_back(e->params);
  }
  return ret;
}

EventLog event_log_from_records(const std::vector<EventRecord>& records) {
  EventLog ret;
  ret.time.reserve(records.size());
  ret.type.reserve(records.size());
  ret.target.reserve(records.size());
  ret.target_index.reserve(records.size());
  ret.requested.reserve(records.size());
  ret.applied.reserve(records.size());
  for (std::vector<EventRecord>::const_iterator r = records.begin();
       r != records.end(); ++r) {
    ret.time.push_back(r->time);
    ret.type.push_back(event_type_to_string(r->type));
    ret.target.push_back(event_target_to_string(r->target));
    ret.target_index.push_back(r->target_index + 1);
    ret.requested.push_back(r->requested);
    ret.applied.push_back(r->applied);
  }
  return ret;
}

}
