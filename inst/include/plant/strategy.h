// -*-c++-*-
#ifndef PLANT_PLANT_STRATEGY_H_
#define PLANT_PLANT_STRATEGY_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <plant/control.h>
#include <plant/internals.h>
#include <plant/util.h>
#include <plant/uniroot.h>
#include <RcppCommon.h> // NA_REAL
#include <plant/uniroot.h>
#include <plant/extrinsic_drivers.h>


namespace plant {

// HEIGHT_INDEX, MORTALITY_INDEX and FECUNDITY_INDEX are a claim about the first
// three state slots of EVERY model, made by around fifty readers across this
// package -- the node, the individual, the stochastic node, and the models
// themselves. Nothing checked it. A model that names those three in another
// order has every one of those readers reach its neighbour's slot, and the
// numbers that come back are finite and plausible: a height that is really a
// mortality is a small positive number.
//
// Called from every model's refresh_indices(), which is where the map that would
// falsify it is built, so it costs one pass per strategy and nothing per
// evaluation.
inline void check_state_layout(const std::map<std::string, int>& state_index,
                               size_t claimed_size, const char* model) {
  // state_size() is a literal on every model and the map is built from
  // state_names(), so the two can disagree. They disagree silently: Internals is
  // sized by the literal, and a name the list has past that length resolves to
  // an index no vector holds.
  if (state_index.size() != claimed_size) {
    util::stop(std::string(model) + " names " +
               util::to_string(static_cast<int>(state_index.size())) +
               " states against the " +
               util::to_string(static_cast<int>(claimed_size)) +
               " its state_size() claims");
  }
  const std::pair<const char*, int> claimed[] = {{"height", HEIGHT_INDEX},
                                                 {"mortality", MORTALITY_INDEX},
                                                 {"fecundity", FECUNDITY_INDEX}};
  for (const std::pair<const char*, int>& c : claimed) {
    const std::map<std::string, int>::const_iterator at =
        state_index.find(c.first);
    if (at == state_index.end()) {
      util::stop(std::string(model) + " names no `" + c.first +
                 "` state, which every generic reader of its state addresses by "
                 "position " + util::to_string(c.second));
    }
    if (at->second != c.second) {
      util::stop(std::string(model) + " names `" + c.first + "` state " +
                 util::to_string(at->second) + ", against the " +
                 util::to_string(c.second) +
                 " every generic reader of it assumes");
    }
  }
}

template <typename E> 
class Strategy {
public:
  typedef E             environment_type;
  typedef std::shared_ptr<Strategy> ptr;

  // update this when the length of state_names changes
  static size_t state_size ();
  // update this when the length of aux_names changes
  size_t aux_size ();

  static std::vector<std::string> state_names();

  // Which of state_names() the model keeps non-negative. The solver rejects a
  // step that lands outside, so the bound holds on the state itself rather than
  // on each reader of it. A strategy with no such state declares none and the
  // check compiles to nothing.
  static std::vector<std::string> non_negative_states() { return {}; }

  std::vector<std::string> aux_names();

  // TODO(#483) : expose this so can access state_names directly
  // In previous attempt couldn't get it to run
  // static std::vector<std::string> state_names() { return strategy_type::state_names(); }
  // the index of variables in the internals extra vector
  std::map<std::string, int> state_index; 
  std::map<std::string, int> aux_index;

  // birth rate spline control points for each species
  // default is constant birth_rate of 1.0
  std::vector<double> birth_rate_x;
  std::vector<double> birth_rate_y = {1.0};
  // whether the spline for each species should be constant fn or not (extrapolation on/off)
  bool is_variable_birth_rate = false;

  bool collect_all_auxiliary;

  void refresh_indices();

  double competition_effect(double size) const;

  double competition_effect_state(Internals<double>& vars);

  void compute_rates(const environment_type& environment, Internals<double>& vars);

  void update_dependent_aux(const int index, Internals<double>& vars);

  // Seed strategy-specific initial ODE states for a newly introduced individual,
  // given its birth environment (called once from Node::compute_initial_conditions
  // before the first compute_rates). Default no-op; strategies that carry an
  // acclimating/tracked state (e.g. TF24f, #525) override this to initialise it at
  // its optimum so there is no birth transient. Resolved on the concrete strategy
  // type by Individual<T,E>, so overriding it here is not required to be virtual.
  void set_initial_states(const environment_type& environment, Internals<double>& vars) {
    (void)environment;
    (void)vars;
  }

  double net_mass_production_dt(const environment_type& environment,
                                double size, double competition_effect_);

  double establishment_probability(const environment_type& environment);

  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  double mortality_dt(double productivity_area, double cumulative_mortality) const;

  double compute_competition(double z, double size) const;

  double initial_size(void) const;

  double size_0;

  // Every Strategy needs a set of Control objects -- these govern
  // things to do with how numerical calculations are performed,
  // rather than the biological control that this class has.
  Control control;

  std::string name;

  ExtrinsicDrivers extrinsic_drivers;
};


}

#endif
