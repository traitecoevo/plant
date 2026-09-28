#include <plant/node_schedule.h>
#include <plant/parameters.h>
#include <plant/util.h>
#include <Rcpp.h>
#include <algorithm> // find, lower_bound, remove
#include <cmath> // log2, exp2
#include <limits> // std::numeric_limits

namespace plant {

NodeSchedule::NodeSchedule(size_t n_species_)
  : n_species(n_species_),
    at(0),
    max_time(std::numeric_limits<double>::infinity()) {
}

size_t NodeSchedule::size() const {
  return schedule.size();
}

size_t NodeSchedule::get_n_species() const {
  return n_species;
}

NodeSchedule NodeSchedule::expand(size_t n_extra,
                                  std::vector<double> times) {
  NodeSchedule ret = *this;
  ret.n_species += n_extra;
  for (size_t i = n_species; i < ret.n_species; ++i) {
    ret.set_times(times, i);
  }
  return ret;
}

void NodeSchedule::clear_times(size_t species_index) {
  std::vector<schedule_entry>::iterator it = schedule.begin();
  while (it != schedule.end()) {
    std::vector<size_t>& species = it->species;
    species.erase(std::remove(species.begin(), species.end(), species_index),
                  species.end());
    // A time no species is introduced at is not a time -- unless something else
    // happens there. Every non-introduction event acts on the patch or the
    // environment, so clearing one species' schedule must not take them with it.
    it = (species.empty() && it->actions.empty()) ? schedule.erase(it) : it + 1;
  }
  reset();
}

double NodeSchedule::get_max_time() const {
  return max_time;
}

std::vector<std::vector<double> > NodeSchedule::get_times() const {
  std::vector<std::vector<double> > ret;
  for (size_t i = 0; i < n_species; ++i) {
    ret.push_back(times(i));
  }
  return ret;
}

void NodeSchedule::set_times(const std::vector<double>& times_,
                               size_t species_index) {
  clear_times(species_index);
  for (std::vector<double>::const_iterator t = times_.begin();
       t != times_.end(); ++t) {
    insert(*t, species_index);
  }
  reset();
}

void NodeSchedule::set_times(const std::vector<std::vector<double> >& times_) {
  util::check_length(times_.size(), n_species);
  for (size_t i = 0; i < n_species; ++i) {
    set_times(times_[i], i);
  }
}

std::vector<double> NodeSchedule::times(size_t species_index) const {
  std::vector<double> ret;
  for (const schedule_entry& i : schedule) {
    if (std::find(i.species.begin(), i.species.end(), species_index) !=
        i.species.end()) {
      ret.push_back(i.time);
    }
  }
  return ret;
}

// The queue flattened back to the wire format: an instant's actions in their
// application order, then its introductions, one Event each.
std::vector<NodeScheduleEvent> NodeSchedule::get_events() const {
  std::vector<Event> ret;
  for (const schedule_entry& e : schedule) {
    for (const Event& a : e.actions) {
      ret.push_back(a);
    }
    for (size_t sp : e.species) {
      ret.push_back(Event(e.time, sp));
    }
  }
  return ret;
}

// Insert one at a time rather than sorting a copy, so the caller's order is
// irrelevant and the grouping and tie-breaking rules live in one place.
// Introductions join their instant's species list; everything else joins its
// actions, held in event_type_rank order and stable within a rank -- two pulses
// at one instant are capped in sequence against the same pool, so the order they
// were given in is the order they must be applied in.
void NodeSchedule::set_all_events(const std::vector<Event>& events_) {
  schedule.clear();
  for (std::vector<Event>::const_iterator e = events_.begin();
       e != events_.end(); ++e) {
    if (e->is_node_introduction()) {
      insert(e->time_introduction(), e->target_index);
    } else {
      std::vector<schedule_entry>::iterator it = entry_at(e->time_introduction());
      std::vector<Event>::iterator a = it->actions.begin();
      while (a != it->actions.end() &&
             event_type_rank(a->type) <= event_type_rank(e->type)) {
        ++a;
      }
      it->actions.insert(a, Event(e->time_introduction(), e->target_index,
                                  e->type, e->target, e->params));
    }
  }
  reset();
}

// Find the instant, or make it. Every entry seeds `times` with its own time, so
// time_introduction() and time_end() are answerable the moment it exists.
std::vector<schedule_entry>::iterator NodeSchedule::entry_at(double time) {
  std::vector<schedule_entry>::iterator it = schedule.begin();
  while (it != schedule.end() && it->time < time) {
    ++it;
  }
  if (it == schedule.end() || !util::identical(it->time, time)) {
    it = schedule.insert(it, schedule_entry{time, {}, {}, {time}});
  }
  return it;
}

void NodeSchedule::reset() {
  at = 0;
}

const schedule_entry& NodeSchedule::next() const {
  if (at >= schedule.size()) {
    Rcpp::stop("All introductions completed");
  }
  return schedule[at];
}

double NodeSchedule::time_end() const {
  if (at >= schedule.size()) {
    Rcpp::stop("All introductions completed");
  }
  return at + 1 < schedule.size() ? schedule[at + 1].time : max_time;
}

void NodeSchedule::pop() {
  if (at >= schedule.size()) {
    Rcpp::stop("Attempt to pop a completed schedule");
  }
  ++at;
}

size_t NodeSchedule::remaining() const {
  return schedule.size() - at;
}

std::vector<odelia::ode::instruction>
NodeSchedule::program_within(double start, double end) const {
  std::vector<odelia::ode::instruction> ret;
  for (const odelia::ode::instruction& s : ode_steps) {
    if (s.time <= start || s.time >= end) {
      continue;
    }
    if (ret.empty()) {
      // The state the replay starts from, which no step reached.
      ret.push_back({start, std::numeric_limits<double>::quiet_NaN()});
    }
    ret.push_back(s);
  }
  return ret;
}

// * R interface
void NodeSchedule::r_clear_times(util::index species_index) {
  clear_times(species_index.check_bounds(n_species));
}

void NodeSchedule::r_set_times(std::vector<double> times_,
                                 util::index species_index) {
  if (!util::is_sorted(times_.begin(), times_.end())) {
    Rcpp::stop("Times must be sorted (increasing)");
  }
  if (times_.size() == 0) {
    Rcpp::stop("Need at least one time");
  }
  if (times_.front() < 0) {
    Rcpp::stop("First time must nonnegative");
  }
  if (times_.back() > max_time) {
    Rcpp::stop("Times cannot be greater than max_time");
  }
  set_times(times_, species_index.check_bounds(n_species));
}

std::vector<double> NodeSchedule::r_times(util::index species_index) const {
  return times(species_index.check_bounds(n_species));
}

void NodeSchedule::r_set_max_time(double x) {
  if (x < 0) {
    Rcpp::stop("max_time must be nonnegative");
  }
  // Guarded on emptiness: this is called on a freshly built schedule, before any
  // time is set, where there is no final scheduled time to be at least.
  if (!schedule.empty() && x < schedule.back().time) {
    Rcpp::stop("max_time must be at least the final scheduled time");
  }
  max_time = x;
  reset();
}

std::vector<double> NodeSchedule::r_ode_times() const {
  std::vector<double> ret;
  ret.reserve(ode_steps.size());
  for (const odelia::ode::instruction& r : ode_steps) {
    ret.push_back(r.time);
  }
  return ret;
}

std::vector<double> NodeSchedule::r_ode_step_sizes() const {
  std::vector<double> ret;
  ret.reserve(ode_steps.size());
  for (const odelia::ode::instruction& r : ode_steps) {
    ret.push_back(r.step_size);
  }
  return ret;
}

// The ODE schedule this run takes: times to stop at, and the sizes that reached
// them where a recorded run is being replayed. Both at once, because a size
// belongs to the time it was recorded with and there is no version of this that
// takes one and keeps the other.
//
// Sizes may be omitted, and that is a grid a caller chose rather than a run a
// caller recorded: the solver steps TO each time instead of by a recorded size.
// A schedule holding either replays it -- there is no flag to turn that on, so to
// integrate freely instead, hand over no schedule.
void NodeSchedule::r_set_ode_steps(std::vector<double> times,
                                   std::vector<double> sizes) {
  if (times.empty()) {
    r_clear_ode_steps();
    return;
  }
  if (times.size() < 2) {
    Rcpp::stop("Need at least two times");
  }
  if (sizes.empty()) {
    // A grid: every time is stepped to, so no size is known anywhere.
    sizes.assign(times.size(), std::numeric_limits<double>::quiet_NaN());
  }
  if (sizes.size() != times.size()) {
    Rcpp::stop("ode_step_sizes must be the same length as ode_times");
  }
  if (!util::identical(times.front(), 0.0)) {
    Rcpp::stop("First time must be exactly zero");
  }
  if (util::is_finite(max_time) && !util::identical(times.back(), max_time)) {
    Rcpp::stop("Last time must be exactly max_time");
  }
  if (!util::is_sorted(times.begin(), times.end())) {
    Rcpp::stop("ode_times must be sorted");
  }
  if (!std::isnan(sizes.front())) {
    Rcpp::stop("First step size must be NaN, the recorded start that no step "
               "reached");
  }
  ode_steps.clear();
  ode_steps.reserve(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    ode_steps.push_back({times[i], sizes[i]});
  }
  if (!util::is_finite(max_time)) {
    max_time = ode_steps.back().time;
  }
  reset();
}

void NodeSchedule::r_clear_ode_steps() {
  ode_steps.clear();
  reset();
}


SEXP NodeSchedule::r_next_introduction() const {
  const schedule_entry& i = next();
  std::vector<size_t> species;
  species.reserve(i.species.size());
  for (const size_t s : i.species) {
    species.push_back(s + 1);
  }
  return Rcpp::List::create(Rcpp::_["time"] = i.time,
                            Rcpp::_["species"] = species,
                            Rcpp::_["time_end"] = time_end());
}

SEXP NodeSchedule::r_all_times() const {
  return Rcpp::wrap(get_times());
}

void NodeSchedule::r_set_all_times(SEXP rx) {
  Rcpp::List x(Rcpp::as<Rcpp::List>(rx));
  // Ensure that we can get all the times out:
  std::vector< std::vector<double> > new_times;
  for (Rcpp::List::iterator el = x.begin(); el != x.end(); ++el) {
    new_times.push_back(Rcpp::as<std::vector<double> >(*el));
  }
  util::check_length(new_times.size(), n_species);
  for (size_t i = 0; i < n_species; ++i) {
    set_times(new_times[i], i);
  }
}

NodeSchedule NodeSchedule::r_copy() const {
  return *this;
}

// * Private methods

// Grouped on exact equality, which is the comparison the run makes when it checks
// it has arrived: a run steps TO a scheduled time and refuses a time that is not
// the one it expected, so two species share an introduction when they share that
// value and not otherwise.
void NodeSchedule::insert(double time, size_t species_index) {
  std::vector<schedule_entry>::iterator it = entry_at(time);
  std::vector<size_t>& species = it->species;
  std::vector<size_t>::iterator s =
    std::lower_bound(species.begin(), species.end(), species_index);
  if (s != species.end() && *s == species_index) {
    // The birth-date coordinate integrates over introduction times and refuses a
    // tie, so this is refused here where it can name what it is rather than
    // later as a duplicate birth date.
    Rcpp::stop("A species cannot be introduced twice at one time");
  }
  species.insert(s, species_index);
}

}
