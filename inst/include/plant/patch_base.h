// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_BASE_H_
#define PLANT_PLANT_PATCH_BASE_H_

#include <numeric>
#include <vector>
#include <odelia/ode_interface.hpp>
#include <plant/control.h>
#include <plant/parameters.h>
#include <plant/util.h>
#include <plant/with_slope.h>

namespace plant {

// What the two patches share:
//   * Patch<T,E>            -- the deterministic patch, whose species carry a
//                              size-density distribution in cohorts; and
//   * StochasticPatch<T,E>  -- the individual-based patch, whose species carry
//                              a finite population of individuals.
//
// Both hold the same parts (parameters, area, environment, species, control)
// and present the same ODE system to the solver: every species' state in
// species order, then the environment's. This CRTP base owns that layout, so it
// is written once: the size, the state, the rates, and the reads and writes of
// them. Two copies of it are what let the stochastic patch drop the
// environment's state for years while the deterministic one carried it, which
// for TF24 meant stochastic runs at frozen soil water.
//
// It is written over the species container, which is the one thing the two
// differ in here, as SpeciesBase is written over the element. What is not here
// is what genuinely differs: how the competition field is built (a trapezium
// over cohorts with a closing boundary interval, against a plain sum over
// individuals), the rates (survival-weighted with an inflow condition, against
// per individual with discrete deaths), introduction, and everything the
// deterministic patch alone does (disturbance, recorded fields, rebinding to
// another scalar, auxiliary variables). As in SpeciesBase, nothing dispatches at
// run time: the hooks below resolve through the CRTP downcast, and a derived
// class that does not declare one gets the base's.
//
// The hooks a derived class may declare (the base calls them on the derived
// type, so a derived declaration hides the default here):
//   * after_set_ode_state()  -- what follows loading a state: building the
//                               field the rates will be taken in, at least.
//   * before_ode_rates()     -- what the rates need computing first, if the
//                               derived class does not compute them on load.
//   * environment_is_given() -- true where the environment is read from
//                               outside, so nothing the species do moves it
//                               and its rates are zero.
//   * check_loaded_nodes()   -- a check on the nodes load_state() introduced,
//                               before their state is loaded.
template <typename Derived, typename T, typename E, typename S>
class PatchBase {
public:
  using value_type = typename T::value_type;

  typedef T               strategy_type;
  typedef E               environment_type;
  typedef Individual<T,E> individual_type;
  typedef S               species_type;
  typedef Parameters<T,E> parameters_type;

  size_t size() const { return species.size(); }
  double time() const { return environment.time; }
  double get_area() const { return area; }

  // Whether any species holds a node: a cohort, or a living individual.
  bool has_nodes() const {
    for (const species_type& s : species) {
      if (s.size() > 0) {
        return true;
      }
    }
    return false;
  }

  const species_type& at_species(size_t species_index) const {
    return species[species_index];
  }

  value_type height_max() const {
    value_type ret = 0.0;
    for (size_t i = 0; i < species.size(); ++i) {
      const value_type h = species[i].height_max();
      if (h > ret) {
        ret = h;
      }
    }
    return ret;
  }

  // The competition profile and its vertical derivative at z, from one pass
  // over the species, each per unit area.
  with_slope<value_type> compute_competition_and_slope(double z) const {
    with_slope<value_type> sum{0.0, 0.0};
    for (size_t i = 0; i < species.size(); ++i) {
      const with_slope<value_type> fs =
        species[i].compute_competition_and_slope(z);
      sum.value += fs.value / area;
      sum.slope += fs.slope / area;
    }
    return sum;
  }
  // The pair's first entry, which each species' pair gives bit for bit. The
  // field is built from the pair (or, deterministically, from its split), never
  // from this, so it is an accessor, and an accessor paying one extra crown
  // evaluation per node is better than a second loop that can disagree with
  // the one above.
  value_type compute_competition(double height) const {
    return compute_competition_and_slope(height).value;
  }

  // Add one species per strategy, each taking this patch's control.
  void add_strategies(std::vector<strategy_type> strategies) {
    for (auto& s : strategies) {
      s.control = control;
      species.push_back(species_type(s));
    }
  }

  // * ODE interface
  size_t ode_size() const {
    return odelia::ode::ode_size(species.begin(), species.end()) +
      environment.ode_size();
  }
  // The nodes' share of the system: zero on a patch with no nodes, whatever
  // state the environment carries.
  size_t node_ode_size() const { return ode_size() - environment.ode_size(); }
  double ode_time() const { return time(); }

  template <typename It> It set_ode_state(It it, double time) {
    it = write_ode_state(it, time);
    d().after_set_ode_state();
    return it;
  }
  template <typename It> It ode_state(It it) const {
    it = odelia::ode::ode_state(species.begin(), species.end(), it);
    return environment.ode_state(it);
  }
  template <typename It> It ode_rates(It it) {
    d().before_ode_rates();
    it = odelia::ode::ode_rates(species.begin(), species.end(), it);
    if (d().environment_is_given()) {
      for (size_t k = 0; k < environment.ode_size(); ++k) {
        *it = 0.0;
        ++it;
      }
      return it;
    }
    return environment.ode_rates(it);
  }

  // * R interface
  parameters_type r_parameters() const { return parameters; }
  environment_type r_environment() const { return environment; }
  std::vector<species_type> r_species() const { return species; }
  species_type r_at(util::index species_index) const {
    return species[species_index.check_bounds(size())];
  }

protected:
  PatchBase(parameters_type p, environment_type e, Control c)
    : parameters(p), area(p.patch_area), environment(e), control(c) {}

  // The state `it` points to, read in the system's layout, and the time; no
  // field is built and no rate computed. set_ode_state() is this and then the
  // derived class's hook.
  template <typename It> It write_ode_state(It it, double time) {
    it = odelia::ode::set_ode_state(species.begin(), species.end(), it);
    it = environment.set_ode_state(it);
    environment.time = time;
    return it;
  }

  // n[i] new nodes in species i, each a copy of the species' next one.
  void introduce_nodes_by_count(const std::vector<size_t>& n) {
    for (size_t i = 0; i < species.size(); ++i) {
      for (size_t j = 0; j < n[i]; ++j) {
        species[i].introduce_new_node();
      }
    }
  }

  // Resets the patch, then puts it at `time` holding n[i] nodes of species i
  // and the flat ODE `state`. Every node is a copy of its species' next one, so
  // what is not in the state (a cohort's birth date) is the same for all.
  void load_state(double time, const std::vector<double>& state,
                  const std::vector<size_t>& n) {
    util::check_length(n.size(), species.size());
    d().reset();
    introduce_nodes_by_count(n);
    util::check_length(state.size(), ode_size());
    d().check_loaded_nodes();
    set_ode_state(state.begin(), time);
  }

  // Each species' state for R, in species order.
  Rcpp::List r_species_states() const {
    Rcpp::List ret;
    for (size_t i = 0; i < species.size(); ++i) {
      ret.push_back(species[i].r_get_state());
    }
    return ret;
  }

  // The uptake of each resource the environment supplies, summed over the
  // species and per unit area, appended to `out`: what closes the
  // environment's own balance (for TF24, soil water). Sized by the resources
  // the environment supplies, which is not its ODE width: TF24 has five soil
  // layers plus four cumulative-flux states that nothing draws on.
  template <typename V>
  void append_resource_depletion(std::vector<V>& out) const {
    for (size_t i = 0; i < environment.n_resources(); i++) {
      const V consumed = std::accumulate(
        species.begin(), species.end(), V(0.0),
        [i](const V& r, const species_type& s) -> V {
          return r + s.consumption_rate(i);
        });
      out.push_back(consumed / area);
    }
  }

  // The defaults for the hooks described above.
  void after_set_ode_state() {}
  void before_ode_rates() {}
  bool environment_is_given() const { return false; }
  void check_loaded_nodes() const {}

  parameters_type parameters;
  double area;
  environment_type environment;
  std::vector<species_type> species;
  Control control;

private:
  Derived&       d()       { return static_cast<Derived&>(*this); }
  const Derived& d() const { return static_cast<const Derived&>(*this); }
};

} // namespace plant

#endif
