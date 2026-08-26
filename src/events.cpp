#include <plant/events.h>
#include <Rcpp.h>

namespace plant {

namespace {
// One table, read in both directions, so a name and its tag cannot drift
// apart. Order follows EventType's, which is also the within-time application
// order.
struct TypeName {
  EventType type;
  const char* name;
  size_t n_params;
};

const TypeName type_names[] = {
  {EventType::RainfallPulse,      "rainfall_pulse",      1},
  {EventType::TemperatureExtreme, "temperature_extreme", 4},
  {EventType::Harvest,            "harvest",             2},
  {EventType::PartialDisturbance, "partial_disturbance", 1},
  {EventType::NodeIntroduction,   "node_introduction",   0}
};

const size_t n_type_names = sizeof(type_names) / sizeof(type_names[0]);
}

EventType event_type_from_string(const std::string& name) {
  for (size_t i = 0; i < n_type_names; ++i) {
    if (name == type_names[i].name) {
      return type_names[i].type;
    }
  }
  std::string known;
  for (size_t i = 0; i < n_type_names; ++i) {
    known += std::string(i > 0 ? ", " : "") + type_names[i].name;
  }
  util::stop("Unknown event type '" + name + "'. Known types: " + known);
  return EventType::NodeIntroduction; // not reached
}

std::string event_type_to_string(EventType type) {
  for (size_t i = 0; i < n_type_names; ++i) {
    if (type == type_names[i].type) {
      return type_names[i].name;
    }
  }
  util::stop("Unknown event type");
  return ""; // not reached
}

size_t event_type_n_params(EventType type) {
  for (size_t i = 0; i < n_type_names; ++i) {
    if (type == type_names[i].type) {
      return type_names[i].n_params;
    }
  }
  util::stop("Unknown event type");
  return 0; // not reached
}

void Events::validate() {
  const size_t n = time.size();
  if (type.size() != n || species_index.size() != n || params.size() != n) {
    util::stop("Events columns must all have the same length (time has " +
               util::to_string(n) + ")");
  }
  for (size_t i = 0; i < n; ++i) {
    if (!util::is_finite(time[i]) || time[i] < 0.0) {
      util::stop("Event " + util::to_string(i + 1) +
                 " has a non-finite or negative time");
    }
    // Also rejects an unknown type name, with the list of known ones.
    const EventType t = event_type_from_string(type[i]);
    const size_t n_expected = event_type_n_params(t);
    if (params[i].size() != n_expected) {
      util::stop("Event " + util::to_string(i + 1) + " (" + type[i] +
                 ") expects " + util::to_string(n_expected) +
                 " parameters but has " + util::to_string(params[i].size()));
    }
  }
}

std::vector<NodeScheduleEvent> to_schedule_events(const Events& events,
                                                  size_t n_species) {
  std::vector<NodeScheduleEvent> ret;
  ret.reserve(events.size());
  for (size_t i = 0; i < events.size(); ++i) {
    const EventType type = event_type_from_string(events.type[i]);
    size_t species = 0;
    if (type == EventType::NodeIntroduction) {
      // 1-based on the way in, as everywhere else in the R interface.
      const size_t raw = events.species_index[i];
      if (raw < 1 || raw > n_species) {
        util::stop("Event " + util::to_string(i + 1) +
                   " has species_index " + util::to_string(raw) +
                   ", outside 1.." + util::to_string(n_species));
      }
      species = raw - 1;
    }
    ret.push_back(NodeScheduleEvent(events.time[i], species, type,
                                    events.params[i]));
  }
  return ret;
}

Events events_from_schedule_events(
    const std::vector<NodeScheduleEvent>& events) {
  Events ret;
  ret.time.reserve(events.size());
  ret.type.reserve(events.size());
  ret.species_index.reserve(events.size());
  ret.params.reserve(events.size());
  for (std::vector<NodeScheduleEvent>::const_iterator e = events.begin();
       e != events.end(); ++e) {
    ret.time.push_back(e->time_introduction());
    ret.type.push_back(event_type_to_string(e->type));
    ret.species_index.push_back(e->species_index + 1); // back to 1-based
    ret.params.push_back(e->params);
  }
  return ret;
}

}
