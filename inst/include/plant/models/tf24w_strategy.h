// -*-c++-*-
#ifndef PLANT_PLANT_TF24W_STRATEGY_H_
#define PLANT_PLANT_TF24W_STRATEGY_H_

#include <plant/models/tf24_strategy.h>

namespace plant {

// TF24w ("w" for water): TF24 in which every cohort carries its own root-zone
// soil column (#681). TF24 integrates one soil column for the patch, which
// every plant reads and drains, so water a plant does not transpire is water
// every plant gets -- and whether a plant that holds water back can ever be
// favoured is a question that model cannot ask. Here the ground under a crown
// is the plant's own: its column takes the rain that falls on it, loses the
// plant's own uptake and its own deep drainage, and exchanges with the rest of
// the patch at the environment's water_mixing_rate. Ground under no crown is
// the environment's existing column. Land moves between columns as crowns
// grow, die and overlap, carrying its water (Patch::private_water_exchange),
// so the patch's water is conserved; at a large mixing rate every column holds
// the patch mean and TF24's water balance is recovered.
//
// Adds one ODE state per soil layer, `water_1` .. `water_N`, volumetric water
// content of the column's layers in the environment's units, appended after
// TF24's six so their indices are unchanged. N is fixed at compile time
// because state_size() is static on the Strategy contract; a patch with more
// layers is refused, one with fewer leaves the spare slots at rate zero.
//
// The physics of a column is the environment's own, called through the const
// helpers TF24_Environment factored out for this (rainfall_floored,
// infiltration_from, soil_K_from_soil_theta, soil_layer_rate,
// psi_from_soil_moist, pulse_capacity). The leaf reads the column's own
// potentials through the virtual soil_water_potential_for_leaf hook.
//
// Limits, enforced where they are first seen: the birth-date size-density
// coordinate (Patch::check_retirement_coordinate), no gradient path
// (SCM::require_birth_date_coordinate) and no stochastic patch
// (StochasticPatch::compute_rates). The regime switches below -- cover above
// or below one -- are `if`s on passive values, which the no-gradient rule
// allows because nothing differentiates through them.
template <typename S = double>
class TF24w_Strategy : public TF24_Strategy<S> {
public:
  using value_type = S;
  using base_type = TF24_Strategy<S>;

  typedef std::shared_ptr<TF24w_Strategy> ptr;
  TF24w_Strategy();

  // Compound scientific version, as TF24f's: "<TF24 version>.<revision>", so a
  // TF24 science change invalidates TF24w too, and this revision moves for
  // changes to the private-water model itself.
  static constexpr int private_water_revision = 1;

  // Soil layers a column carries: the environment's default. The name is what
  // the HasPrivateWater concept tests for.
  static constexpr size_t n_water_layers = 5;

  static size_t state_size() { return base_type::state_size() + n_water_layers; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = base_type::state_names();
    for (size_t i = 0; i < n_water_layers; ++i) {
      ret.push_back("water_" + std::to_string(i + 1));
    }
    return ret;
  }
  // A column's layers, like the storage pool, hold exact flows inside their
  // range: a negative content is a step that overshot, rejected and retried.
  static std::vector<std::string> non_negative_states() {
    std::vector<std::string> ret = base_type::non_negative_states();
    for (size_t i = 0; i < n_water_layers; ++i) {
      ret.push_back("water_" + std::to_string(i + 1));
    }
    return ret;
  }
  // Three readings beside TF24's: the column's infiltration and deep drainage
  // (m yr^-1 per m2 of its footprint), which the exchange pass area-weights
  // into the environment's accumulators, and the crown area (m2).
  std::vector<std::string> aux_names() override {
    std::vector<std::string> ret = base_type::aux_names();
    ret.push_back("water_infiltration");
    ret.push_back("water_drainage");
    ret.push_back("crown_area");
    return ret;
  }

  void refresh_indices();
  void compute_rates(const TF24_Environment<S>& environment, Internals<S>& vars);
  // Chains to TF24's, then seeds the column at the content of the ground a
  // newborn takes over (the environment caches it): the bare ground while
  // crowns leave some, the private columns' mean once they cover it. That is
  // what keeps an introduction from making or losing water.
  void set_initial_states(const TF24_Environment<S>& environment, Internals<S>& vars);
  // The no-state overload reads the same content.
  S establishment_probability(const TF24_Environment<S>& environment);
  using base_type::establishment_probability;

  // The ground a plant's roots own: the area its crown covers, with L_crown
  // the leaf area a crown carries per unit of that ground.
  S crown_area(const S& area_leaf_) const { return area_leaf_ / this->pars.L_crown; }

  // Cached slots, resolved in refresh_indices().
  int state_idx_water_first = -1;
  int aux_idx_water_infiltration = -1;
  int aux_idx_water_drainage = -1;
  int aux_idx_crown_area = -1;

protected:
  const std::vector<S>& soil_water_potential_for_leaf(
      const TF24_Environment<S>& environment) override {
    (void)environment;
    return psi_soil_own_;
  }
  // The column's potentials for the leaf, filled before the base rates run.
  std::vector<S> psi_soil_own_;
  std::vector<S> water_flux_own_;
};

template <typename S>
typename TF24w_Strategy<S>::ptr make_strategy_ptr(TF24w_Strategy<S> s);

template <typename S>
TF24w_Strategy<S>::TF24w_Strategy() {
  this->name = "TF24w";
  this->refresh_indices();
}

template <typename S>
void TF24w_Strategy<S>::refresh_indices() {
  // Rebuilds both maps from the virtual aux_names(), so the three readings
  // above are registered too, and caches TF24's own slots.
  base_type::refresh_indices();
  state_idx_water_first = static_cast<int>(base_type::state_size());
  for (size_t i = 0; i < n_water_layers; ++i) {
    this->state_index["water_" + std::to_string(i + 1)] =
      state_idx_water_first + static_cast<int>(i);
  }
  aux_idx_water_infiltration = this->aux_index.at("water_infiltration");
  aux_idx_water_drainage = this->aux_index.at("water_drainage");
  aux_idx_crown_area = this->aux_index.at("crown_area");
  check_state_layout(this->state_index, state_size(), "TF24w");
}

template <typename S>
void TF24w_Strategy<S>::compute_rates(const TF24_Environment<S>& environment,
                                      Internals<S>& vars) {
  using odelia::util::to_passive;
  const int n = environment.get_soil_number_of_depths();
  if (n > static_cast<int>(n_water_layers)) {
    util::stop("TF24w carries " + util::to_string(n_water_layers) +
               " soil layers per cohort and this environment has " +
               util::to_string(n));
  }
  // The leaf solves against this column.
  psi_soil_own_.resize(n);
  for (int i = 0; i < n; ++i) {
    psi_soil_own_[i] =
      environment.psi_from_soil_moist(vars.state(state_idx_water_first + i), i);
  }
  // TF24's rates: the leaf solve, growth, mortality, fecundity, storage, and
  // the per-layer uptake in consumption_rates (m3 yr^-1 per individual). The
  // patch reads that uptake for its accumulator and passes the environment's
  // column none of it.
  base_type::compute_rates(environment, vars);

  const S& area_leaf_ = vars.aux(this->aux_idx_competition_effect);
  const S area = crown_area(area_leaf_);
  vars.set_aux(aux_idx_crown_area, area);
  // Where crowns overlap, each plant's ground is its crown's share of the
  // patch; the cover is the exchange pass's, formed from the state before any
  // column was rated.
  const S area_eff = area / std::max(1.0, environment.private_water_cover());

  const double rainfall = environment.rainfall_floored();
  const S infiltration =
    environment.infiltration_from(rainfall, vars.state(state_idx_water_first));
  const double m = environment.water_mixing_rate;
  water_flux_own_.resize(n);
  S water_input;
  for (int i = 0; i < n; ++i) {
    const int slot = state_idx_water_first + i;
    const S& theta = vars.state(slot);
    water_input = (i == 0) ? infiltration : water_flux_own_[i - 1];
    water_flux_own_[i] = environment.soil_K_from_soil_theta(theta, i);
    // Uptake per unit of the plant's own ground. Crown area is leaf area over
    // L_crown and uptake is per leaf area times leaf area, so this is leaf-area
    // density times uptake per leaf area, the same at every size.
    const S uptake = vars.consumption_rate(i) / area_eff;
    S rate = environment.soil_layer_rate(water_input, water_flux_own_[i], uptake,
                                         theta, i);
    if (m > 0.0) {
      rate += m * (S(environment.private_water_mean(i)) - theta);
    }
    vars.set_rate(slot, rate);
  }
  for (int i = n; i < static_cast<int>(n_water_layers); ++i) {
    vars.set_rate(state_idx_water_first + i, S(0.0));
  }
  vars.set_aux(aux_idx_water_infiltration, infiltration);
  vars.set_aux(aux_idx_water_drainage, water_flux_own_[n - 1]);
}

template <typename S>
void TF24w_Strategy<S>::set_initial_states(const TF24_Environment<S>& environment,
                                           Internals<S>& vars) {
  base_type::set_initial_states(environment, vars);
  const int n = environment.get_soil_number_of_depths();
  for (int i = 0; i < static_cast<int>(n_water_layers); ++i) {
    vars.set_state(state_idx_water_first + i,
                   i < n ? S(environment.private_water_newborn(i)) : S(0.0));
  }
}

template <typename S>
S TF24w_Strategy<S>::establishment_probability(const TF24_Environment<S>& environment) {
  const int n = environment.get_soil_number_of_depths();
  psi_soil_own_.resize(n);
  for (int i = 0; i < n; ++i) {
    psi_soil_own_[i] =
      environment.psi_from_soil_moist(S(environment.private_water_newborn(i)), i);
  }
  return base_type::establishment_probability(environment);
}

template <typename S>
typename TF24w_Strategy<S>::ptr make_strategy_ptr(TF24w_Strategy<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24w_Strategy<S> >(s);
}

}

#endif
