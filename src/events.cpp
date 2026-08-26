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
  bool accepts_species;
};

const TypeInfo type_info[] = {
  // Soil water belongs to the patch, not to any one species, so a pulse
  // cannot be narrowed to one.
  {EventType::RainfallPulse,    "rainfall_pulse",    1,
   EventTarget::Environment, false},
  // temperature, duration, temperature_crit, sensitivity
  {EventType::HeatDamage,       "heat_damage",       4,
   EventTarget::Patch,       true},
  // fraction, height_min, height_max
  {EventType::Thinning,         "thinning",          3,
   EventTarget::Patch,       true},
  {EventType::NodeIntroduction, "node_introduction", 0,
   EventTarget::Species,     true}
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
  return info_for(type).accepts_species;
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
    if (tg == EventTarget::Species && !event_type_accepts_species(t)) {
      util::stop(at + " cannot be aimed at a single species");
    }
    if (tg != EventTarget::Species && t == EventType::NodeIntroduction) {
      util::stop(at + " must name the species being introduced");
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
    size_t index = 0;
    if (target == EventTarget::Species) {
      // 1-based on the way in, as everywhere else in the R interface.
      const size_t raw = events.target_index[i];
      if (raw < 1 || raw > n_species) {
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
