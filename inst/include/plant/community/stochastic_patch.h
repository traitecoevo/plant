// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_PATCH_H_
#define PLANT_PLANT_STOCHASTIC_PATCH_H_

#include <plant/parameters.h>
#include <plant/community/patch_system.h>
#include <plant/solvers/schedule.h>
#include <plant/solvers/events.h>
#include <plant/community/stochastic_species.h>
#include <plant/util.h>
#include <plant/numerics/with_slope.h>

namespace plant {

// NOTE: compute_environment() here might fail (especially for
// rare seed arrivals) because the adaptive refinement can't deal with
// the sharp corners that are implied.  The simplest thing to do is to
// tone down the tolerance (fast_control() seems good enough) but that
// might not be enough.  It might be best to do something more clever
// than just fail on refining, but doing that will require
// confirmation that the issue is simply in a couple of places rather
// than throughout.  Running the spline piecewise would be the best
// bet there.
template <typename T, typename E>
class StochasticPatch
  : public PatchSystem<StochasticPatch<T,E>, T, E, StochasticSpecies<T,E>> {
  typedef PatchSystem<StochasticPatch<T,E>, T, E, StochasticSpecies<T,E>> base_type;
  friend base_type;
public:
  using typename base_type::value_type;
  using typename base_type::strategy_type;
  using typename base_type::environment_type;
  using typename base_type::individual_type;
  using typename base_type::species_type;
  using typename base_type::parameters_type;

  StochasticPatch(parameters_type p, environment_type e, Control c);
  void reset();

  // The patch's time is wherever the solver last evaluated it, which on a
  // leg's final stage can round one ulp away from where the leg landed. The
  // runner sets it to the landing time, so a re-sync does not move the solver.
  void set_time(double t) { environment.time = t; }

  bool introduce_new_node(size_t species_index);
  void introduce_new_node_and_update(size_t species_index);

  std::vector<size_t> deaths();
  // A scheduled event other than a death check or an introduction, which the
  // runner applies itself. The same contract as Patch::apply_event, drawn per
  // individual: a harvest removes each living individual in its size band with
  // probability `fraction` (1 removes them all), a climate extreme each in its
  // target with probability 1 - survival, and a resource pulse acts on the
  // environment as it does there. A removal reports {fraction applied,
  // individuals drawn for, individuals removed per unit area}.
  EventRecord apply_event(const ScheduleEvent& event);

  // * R interface
  void r_set_state(double time,
                   const std::vector<double>& state,
                   const std::vector<size_t>& n) {
    this->load_state(time, state, n);
  }
  Rcpp::List r_get_state() const;
  bool r_introduce_new_node(util::index species_index) {
    return introduce_new_node(species_index.check_bounds(this->size()));
  }
  void r_introduce_new_node_and_update(util::index species_index) {
    introduce_new_node_and_update(species_index.check_bounds(this->size()));
  }

  // These are only here because they wrap private functions.
  void r_compute_environment() {compute_environment();}
  void r_compute_rates() {compute_rates();}
private:
  using base_type::parameters;
  using base_type::area;
  using base_type::environment;
  using base_type::species;
  using base_type::control;

  void compute_environment();
  void compute_rates();

  // PatchSystem's hooks. A loaded state is the rates' state too: they are
  // computed here, where the deterministic patch computes them when they are
  // asked for. The environment is always this patch's own, and an individual
  // carries nothing that loading its state could leave inconsistent.
  void after_set_ode_state() {
    compute_environment();
    compute_rates();
  }
  void before_ode_rates() {}
  bool environment_is_given() const { return false; }
  void check_loaded_nodes() const {}

  // Scratch for compute_rates(), which the solver calls on every derivatives
  // evaluation. A member reserved once in reset() and cleared (not freed) after
  // use, so the RHS does not allocate.
  std::vector<double> resource_depletion;
  // Scratch for compute_environment(): each species' living individuals,
  // tallest first, taken once per profile and read at every height on it.
  std::vector<std::vector<const typename species_type::node_type*>> tallest_first;
};

template <typename T, typename E>
StochasticPatch<T,E>::StochasticPatch(parameters_type p, environment_type e, Control c)
  : base_type(p, e, c, c.stochastic_patch_area) {
  if (!(area > 0.0 && std::isfinite(area))) {
    util::stop("stochastic_patch_area must be positive and finite");
  }
  parameters.validate();
  this->add_strategies(parameters.strategies);
  reset();
}

template <typename T, typename E>
void StochasticPatch<T,E>::reset() {
  for (auto& s : species) {
    s.clear();
  }
  resource_depletion.reserve(environment.n_resources());
  environment.clear();
  compute_environment();
  compute_rates();
}

template <typename T, typename E>
void StochasticPatch<T,E>::compute_environment() {
  if (this->height_max() > 0.0) {
    // Written as std::vector<double> this still compiles, taking the value of an
    // active profile, and the field's knot values and slopes would then be
    // constants with nothing raised to say so.
    //
    // Individuals rather than nodes here, so there is no trapezium over a grid and
    // nothing to accumulate across knots: the per-knot sum IS the reduction, and
    // this walks the knots itself. The sum is the patch's
    // compute_competition_and_slope(), with each species ordered once.
    tallest_first.resize(species.size());
    for (size_t i = 0; i < species.size(); ++i) {
      species[i].order_by_height(tallest_first[i]);
    }
    auto f = [&] (const std::vector<double>& x, std::vector<value_type>& y,
                  std::vector<value_type>& m) -> void {
      for (size_t k = 0; k < x.size(); ++k) {
        with_slope<value_type> sum{0.0, 0.0};
        for (size_t i = 0; i < species.size(); ++i) {
          const with_slope<value_type> fs =
            species[i].compute_competition_and_slope(x[k], tallest_first[i]);
          sum.value += fs.value / area;
          sum.slope += fs.slope / area;
        }
        y[k] = sum.value;
        m[k] = sum.slope;
      }
    };
    environment.compute_environment(f, this->height_max());
    // Kept for its capacity; the pointers would outlive a change to the plants.
    for (auto& order : tallest_first) {
      order.clear();
    }
  } else {
    environment.clear_environment();
  }
}


template <typename T, typename E>
void StochasticPatch<T,E>::compute_rates() {
  if constexpr (HasPrivateWater<T>) {
    util::stop("TF24w's private root-zone water is not implemented for the "
               "stochastic patch (#681)");
  }
  for (size_t i = 0; i < this->size(); ++i) {
    species[i].compute_rates(environment);
  }
  this->append_resource_depletion(resource_depletion);
  environment.compute_rates(resource_depletion);
  resource_depletion.clear();
}

// In theory, this could be done more efficiently by, in the introdudce_new_node
// case, using the values stored in the species offspring.  But we don't
// really get that here.  It might be better to move introdudce_new_node /
// introduce_new_node_and_update within Species, given this.
template <typename T, typename E>
void StochasticPatch<T,E>::introduce_new_node_and_update(size_t species_index) {
  // Add a offspring, setting ODE variables based on the *current* light environment
  species[species_index].introduce_new_node(environment);
  // Then we update the light environment.
  compute_environment();
}

template <typename T, typename E>
bool StochasticPatch<T,E>::introduce_new_node(size_t species_index) {
  // The draw decides whether an individual is added, so the probability is read
  // at its value here.
  const double pr_germinate = odelia::util::to_passive(
    species[species_index].establishment_probability(environment));
  const bool added = unif_rand() < pr_germinate;
  if (added) {
    introduce_new_node_and_update(species_index);
  }
  return added;
}

template <typename T, typename E>
std::vector<size_t> StochasticPatch<T,E>::deaths() {
  std::vector<size_t> ret;
  ret.reserve(this->size());
  bool recompute = false;
  for (auto& s : species) {
    const size_t n_deaths = s.deaths();
    ret.push_back(n_deaths);
    recompute = recompute || n_deaths > 0;
  }
  if (recompute) {
    compute_environment();
    compute_rates();
  }
  return ret;
}

template <typename T, typename E>
EventRecord StochasticPatch<T,E>::apply_event(const ScheduleEvent& event) {
  EventRecord rec;
  rec.time = this->time();
  rec.type = event.type;
  rec.target = event.target;
  rec.target_index = event.target_index;
  rec.requested = event.params;

  std::vector<size_t> targets;
  if (event.target == EventTarget::Species) {
    targets.push_back(event.target_index);
  } else {
    for (size_t i = 0; i < species.size(); ++i) {
      targets.push_back(i);
    }
  }
  auto remove = [&](double fraction, auto pr) {
    size_t drawn = 0, removed = 0;
    for (size_t i : targets) {
      const auto r = species[i].remove(pr);
      drawn += r.first;
      removed += r.second;
    }
    rec.applied = {fraction, static_cast<double>(drawn),
                   static_cast<double>(removed) / area};
    compute_environment();
    compute_rates();
  };

  switch (event.type) {
  case EventType::DeathCheck:
  case EventType::NodeIntroduction:
    util::stop("Death checks and introductions are the stochastic runner's "
               "to apply");
    break;

  case EventType::ResourcePulse:
    if (event.target_index >= environment.n_resources()) {
      util::stop("Resource " + util::to_string(event.target_index + 1) +
                 " does not exist: this environment has " +
                 util::to_string(environment.n_resources()) + " resources");
    }
    rec.applied = environment.add_resource_pulse(event.target_index,
                                                 event.params.at(0));
    // The environment's rates read its pools.
    compute_rates();
    break;

  case EventType::Harvest: {
    const double fraction = event.params.at(0);
    const double size_min = event.params.at(1);
    const double size_max = event.params.at(2);
    if (!(fraction >= 0.0 && fraction <= 1.0)) {
      util::stop("Harvest fraction must be in [0, 1]");
    }
    remove(fraction, [&](double size) {
      return (size >= size_min && size <= size_max) ? fraction : 0.0;
    });
    break;
  }

  case EventType::ClimateExtreme: {
    const double duration = event.params.at(1);
    if (!util::is_finite(duration) || duration < 0.0) {
      util::stop("Climate extreme duration must be finite and non-negative");
    }
    const double phi = climate_extreme_survival(
      event.params.at(0), duration, event.params.at(2), event.params.at(3));
    remove(1.0 - phi, [&](double) { return 1.0 - phi; });
    break;
  }
  }
  return rec;
}

template <typename T, typename E>
Rcpp::List StochasticPatch<T, E>::r_get_state() const {
  return Rcpp::List::create(_["time"] = this->time(),
                            _["species"] = this->r_species_states(),
                            _["env"] = environment.r_get_state());
}

}

#endif
