// Built from  inst/include/plant/models/ff16_environment.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_ENVIRONMENT_H_
#define PLANT_PLANT_TF24_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <odelia/interpolator.hpp>
#include <limits>
#include <cmath>

using namespace Rcpp;

namespace plant {

class TF24_Environment : public Environment {
public:
  std::vector<double> resolve_soil_parameter_values(SEXP values,
                                                    int n,
                                                    double default_value,
                                                    const std::string &name) const {
    if (Rf_isNull(values)) {
      return std::vector<double>(n, default_value);
    }

    const auto out = Rcpp::as<std::vector<double> >(values);
    if (out.size() != static_cast<size_t>(n)) {
      throw std::invalid_argument(
        name + " must have length equal to soil_number_of_depths when provided.");
    }

    return out;
  }

  double soil_parameter_value(const std::vector<double> &layer_values,
                              double fallback,
                              size_t layer) const {
    if (use_layered_soil_parameters &&
        layer_values.size() == static_cast<size_t>(soil_number_of_depths)) {
      return layer_values[layer];
    }
    return fallback;
  }

  // constructor for R interface - default settings can be modified
  // except for soil_number_of_depths
  // which are only updated on construction
  
  TF24_Environment(bool light_availability_spline_rescale_usually = true,
                   int soil_number_of_depths = 5, 
                   double delta_z = 9999, // not using this
                   double soil_moist_sat = 0.428, // saturated soil moisture content (m3 water m^-3 soil) 
                   double K_sat = 163.0411, //saturated hydraulic conductivity of soil
                   double a_psi = 1.78e3, // not currently being used
                   double n_psi = 6.57, // not currently being used
                   double a_infil = 1, // infiltration switch (0-1), 0 no runoff, 1 runoff
                   double b_infil = 8, // unitless, determines infiltration rate
                   double depth = 1.5)  // total depth of soil (m)
      : delta_z(delta_z),
      soil_moist_sat(soil_moist_sat),
      a_psi(a_psi),
      n_psi(n_psi),
      K_sat(K_sat),
      a_infil(a_infil),
      b_infil(b_infil),
      depth(depth)
  {
    time = 0.0;



    // Shading defaults have lower tolerance which are overwritten for speed
    light_availability = ResourceSpline(
                   1e-4,  // light_availability_spline_tol,
                   17,    // light_availability_spline_nbase,
                   16,    // light_availability_spline_max_depth,
                   true //light_availability_spline_rescale_usually)
                  );

    ExtrinsicDrivers extrinsic_drivers;

    extrinsic_drivers_set_constant("PPFD",1800);
    extrinsic_drivers_set_constant("rainfall",1);
    extrinsic_drivers_set_constant("atm_vpd",1);
    extrinsic_drivers_set_constant("ca",40);
    extrinsic_drivers_set_constant("leaf_temp",25);
    extrinsic_drivers_set_constant("atm_o2_kpa",21);
    // 101.3 kPa, standard sea-level pressure -- and, more to the point, the value
    // the leaf model's photosynthesis side has always assumed. Its ppm -> Pa
    // conversion was the hard-coded constant 0.1013 = 1e-6 * 101300 Pa, i.e. 101.3
    // kPa in disguise, while this driver said 100.5. The two disagreed, so the
    // conductance side of the model responded to 100.5 while Gamma*, Kc, Ko, Km and
    // the ci root-find bounds silently assumed 101.3.
    //
    // The leaf package now derives the conversion from atm_kpa (phylloptim #15 item
    // 10c), which makes the model self-consistent at any pressure -- and turns the
    // disagreement into a 2.4% shift in TF24 output. **This line is why**: 100.5
    // arrived in `34d46ac2` ("Simplify scm & environment interface", #446), a pure
    // interface refactor that does not mention atmospheric pressure, and no
    // rationale for it is recorded anywhere. Every leaf-level test uses 101.3.
    // So it was an artefact, not a site elevation, and pinning it to the value the
    // rest of the model already assumed keeps the published numbers instead of
    // re-baselining them against an accident. Set it per-site if you mean altitude.
    extrinsic_drivers_set_constant("atm_kpa",101.3);
    extrinsic_drivers_set_constant("wind_speed",2.0); // U0, m s^-1 (Penman-Monteith, #523)

    set_soil_number_of_depths(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, soil_moist_sat*0.5));
  };
  
  // Number of cumulative auxilliary variables to track in soil moisture model
  // Trailing state slots that accumulate diagnostic fluxes rather than
  // participating in the water balance: sum_rainfall, sum_infiltration,
  // sum_drainage, sum_resource_depletion, sum_pulse_runoff.
  //
  // The fifth is fed only by rainfall-pulse events (#522). That is deliberate,
  // not an oversight: the ODE controller takes its error norm over *every*
  // state component, so an accumulator with a non-zero rate would join the
  // step-size decision and could move existing TF24 trajectories. Held at rate
  // zero it contributes exactly zero to the norm, so a run with no pulses is
  // bit-identical to one from before this slot existed. Accumulating the
  // continuous saturation-excess runoff here too would be a deliberate change
  // to what TF24 reports, and needs its own decision.
  static constexpr size_t aux_num = 5;
  
  // Setup soil water distribtuion
  void set_soil_number_of_depths(int n) {
    soil_number_of_depths = n;
    
    vars = Internals(soil_number_of_depths + aux_num);
    initial_states = vars.states;

    z.resize(soil_number_of_depths);
    z_mid.resize(soil_number_of_depths);
    dz.resize(soil_number_of_depths);
    // positive downwards
    water_flux.resize(soil_number_of_depths);

    delta_z = depth / soil_number_of_depths;

    for (int i = 0; i < soil_number_of_depths; i++)
    {
      z[i] = (i + 1) * delta_z;
      if (i == 0) {
        z_mid[i] = z[i] / 2.0;
      } else {
        z_mid[i] = (z[i - 1] + z[i]) / 2.0;
      }
    }

    for (int i = 0; i < soil_number_of_depths; i++)
    {
      dz[i] = delta_z;
    }

    psi_soil_cache_.resize(soil_number_of_depths);
    psi_soil_cache_state_.resize(soil_number_of_depths);
    psi_soil_cache_valid_ = false;

    use_layered_soil_parameters = false;
    soil_moist_sat_layers.clear();
    K_sat_layers.clear();
    a_psi_layers.clear();
    n_psi_layers.clear();
  }

  void set_soil_parameters(int n,
                           SEXP soil_moist_sat_values,
                           SEXP K_sat_values,
                           SEXP a_psi_values,
                           SEXP n_psi_values) {
    set_soil_number_of_depths(n);

    soil_moist_sat_layers = resolve_soil_parameter_values(
      soil_moist_sat_values, n, soil_moist_sat, "soil_moist_sat");
    K_sat_layers = resolve_soil_parameter_values(
      K_sat_values, n, K_sat, "K_sat");
    a_psi_layers = resolve_soil_parameter_values(
      a_psi_values, n, a_psi, "a_psi");
    n_psi_layers = resolve_soil_parameter_values(
      n_psi_values, n, n_psi, "n_psi");

    use_layered_soil_parameters = true;
    psi_soil_cache_valid_ = false;
  }

  int get_soil_number_of_depths() const {return soil_number_of_depths;}

  // Water is taken up per soil layer; the trailing aux_num state slots only
  // accumulate diagnostics and are never consumed.
  size_t n_resources() const override {
    return static_cast<size_t>(soil_number_of_depths);
  }

  std::vector<double> get_soil_mid_depths() const { return z_mid; }

  // TODO: should we use auxilliary in internals
  std::vector<double> water_flux;
  std::vector<double> z;
  std::vector<double> z_mid;
  std::vector<double> dz;
  mutable std::vector<double> psi_soil_cache_;
  mutable std::vector<double> psi_soil_cache_state_;
  mutable bool psi_soil_cache_valid_ = false;

  // Per-driver memo for the time-varying extrinsic drivers (see get_* above).
  // cache_time NaN-initialised so the first call always misses (NaN != time).
  double cached_driver_(const std::string &name, double &cache_val,
                        double &cache_time) const {
    if (cache_time != time) {
      cache_val = extrinsic_drivers.evaluate(name, time);
      cache_time = time;
    }
    return cache_val;
  }
  static constexpr double NAN_TIME_ = std::numeric_limits<double>::quiet_NaN();
  mutable double ppfd_cache_ = 0, atm_vpd_cache_ = 0, ca_cache_ = 0,
                 leaf_temp_cache_ = 0, atm_o2_kpa_cache_ = 0, atm_kpa_cache_ = 0,
                 wind_speed_cache_ = 0;
  mutable double ppfd_cache_time_ = NAN_TIME_, atm_vpd_cache_time_ = NAN_TIME_,
                 ca_cache_time_ = NAN_TIME_, leaf_temp_cache_time_ = NAN_TIME_,
                 atm_o2_kpa_cache_time_ = NAN_TIME_, atm_kpa_cache_time_ = NAN_TIME_,
                 wind_speed_cache_time_ = NAN_TIME_;

  // The soil state a run begins from: whatever set_soil_water_state last set,
  // which the R interface lets a caller choose. Restored by clear_state().
  std::vector<double> initial_states;

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline light_availability;

  // Light interface
  bool canopy_rescale_usually;
  //distance between layers
  int soil_number_of_depths;
  double delta_z;

  double depth;
  //saturated soil moisture
  double soil_moist_sat;
  //Saturated soil hydraulic conductivity
  double K_sat;
  double a_psi;
  double n_psi;
  std::vector<double> soil_moist_sat_layers;
  std::vector<double> K_sat_layers;
  std::vector<double> a_psi_layers;
  std::vector<double> n_psi_layers;
  bool use_layered_soil_parameters = false;
  double a_infil;
  double b_infil;

  // Residual soil moisture floor (m3 m^-3). The water balance cannot dry a
  // layer below this (see the positivity guard in compute_rates), and
  // psi_from_soil_moist is evaluated at this floor so the retention curve
  // (which diverges as theta->0) stays finite. Well below any realistic
  // operating moisture, so it does not perturb non-drought runs. See issue
  // #485.
  double soil_moist_residual = 1e-2;

  // Ceiling on soil water potential (MPa), issue #549. Near the residual floor
  // the retention curve returns psi ~ 1e8 MPa, which drives the leaf hydraulic
  // solve non-finite and blows up the soil feedback. Root conductance has
  // already fallen to ~0 by a few MPa, so clamping here changes nothing
  // physically while keeping the numerics finite under extreme drought.
  double soil_psi_max_ = 1e3;

  // Ability to prescribe a fixed value
  // TODO: add setting to set other variables like water
  void set_fixed_environment(double value, double height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return light_availability.get_value_at_height(height);
  }

  // Highest height covered by the light spline; hoist out of hot per-point
  // loops and feed back into the capped get_environment_at_height() overload.
  double max_environment_height() const {
    return light_availability.max_height();
  }

  double get_environment_at_height(double height, double cap) const {
    return light_availability.get_value_at_height(height, cap);
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }

  // ------------------------------------------------------------------
  // SOIL WATER BALANCE (multi-layer bucket model)
  // ------------------------------------------------------------------
  // Each soil layer i is a bucket holding volumetric moisture vars.state(i)
  // (m3 water m^-3 soil). The rate of change is a simple mass balance:
  //
  //   d(theta_i)/dt = (water_in_i - water_out_i - root_uptake_i) / dz[i]
  //
  // where:
  //   * water_in_0   = infiltration (rainfall reduced by a saturation-excess
  //                    runoff term controlled by a_infil/b_infil);
  //   * water_in_i>0 = drainage out of the layer above (water_flux[i-1]);
  //   * water_out_i  = gravitational drainage = soil_K_from_soil_theta(theta_i),
  //                    a Clapp & Hornberger (1978) / Zeng & Decker (2009)
  //                    unsaturated hydraulic conductivity;
  //   * root_uptake_i= resource_depletion[i], supplied by the plants via the
  //                    strategy's evapotranspiration_dt (m yr^-1).
  //
  // The final `aux_num` state slots accumulate diagnostic cumulative fluxes
  // (rainfall, infiltration, deep drainage, total root uptake).
  //
  // This is an explicit, first-order representation; drainage is instantaneous
  // single-direction (no upward capillary flux between layers - that is handled
  // hydraulically inside the plant via E_from_Soil_to_Root_Collar).
  virtual void compute_rates(std::vector<double> const &resource_depletion)
  {

    double water_input;
    // Rainfall is floored at zero. Drivers are interpolated with a cubic
    // spline, which overshoots badly on intermittent forcing: a realistic daily
    // series with a ~10% wet-day fraction evaluates negative at ~45% of points,
    // reaching -5.7 m yr^-1. Unfloored, negative rainfall gives negative
    // infiltration and a negative layer-0 rate, which is unphysical (rain
    // drying the soil) and fails in two different ways depending on wetness:
    // above residual moisture the water really is removed, while at/below
    // residual the guard further down clamps the rate to zero, so the removal
    // is recorded in `sum_rainfall` but never applied and the water budget
    // stops closing. Drylands sit at residual for much of the year, so that
    // second case is the common one in the intended application.
    //
    // Flooring removes both. It is a bound on the sign, NOT a correction to the
    // interpolation: the spline conserves the integral exactly (undershoot is
    // compensated by overshoot), so discarding the negative lobes raises total
    // rainfall by the undershoot area -- ~7% for the series above. The real
    // remedy is not to spline an intermittent series; run
    // `check_driver_interpolation()` on any such driver, which reports the
    // undershoot area and warns.
    double rainfall = std::max(0.0, extrinsic_drivers.evaluate("rainfall", time));
    const double soil_moist_sat_0 =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, 0);
    double infiltration = rainfall * std::max(
      0.0,
      1 - a_infil * std::pow(vars.state(0) / soil_moist_sat_0, b_infil));
    double total_resource_depletion = 0;


    // treat each soil layer as a separate resource pool
    for (size_t i = 0; i < soil_number_of_depths; i++)
    {

      // initial representation of drainage; to be improved
      if (i == 0)
      {
        water_input = infiltration;
      }
      else
      {
        // m3 m^-2
        water_input = water_flux[i-1];
      }
        // TODO: m3 m^-2
      water_flux[i] = soil_K_from_soil_theta(vars.state(i), i);
      // this function does runoff

      // Positivity guard (issue #485): a layer at or below the residual
      // moisture theta_r is not dried further (only rewetting is allowed). This
      // keeps the explicit fixed-step solver from driving a drought-stressed
      // layer to theta <= 0, where the retention curve psi_from_soil_moist and
      // the conductivity curve soil_K_from_soil_theta go non-finite. Wetter
      // layers are unaffected, so non-drought runs are unchanged.
      const double theta = vars.state(i);
      double rate = (water_input - water_flux[i] - resource_depletion[i]) / dz[i];
      // Positivity guard (issue #485), hardened for #549: at/below the residual
      // moisture a layer must not be dried further. The original `rate < 0.0`
      // test let a *non-finite* rate through, because `NaN < 0.0` is false in
      // IEEE 754 -- a NaN soil_consumption_ (from the retention curve's enormous
      // psi_soil near theta_r) then wrote straight into the soil state. `!(rate
      // > 0.0)` is true for NaN and for rate <= 0, so only genuine rewetting is
      // allowed and NaN/negative rates are clamped to 0.
      if (theta <= soil_moist_residual && !(rate > 0.0)) {
        rate = 0.0;
      }
      vars.set_rate(i, rate);
      total_resource_depletion += resource_depletion[i];
    }
      vars.set_rate(soil_number_of_depths, rainfall);
      vars.set_rate(soil_number_of_depths + 1, infiltration);
      vars.set_rate(soil_number_of_depths + 2, water_flux[soil_number_of_depths - 1]);
      vars.set_rate(soil_number_of_depths + 3, total_resource_depletion);
      // Pulse runoff has no continuous source: it only ever moves in
      // add_water_pulse(). Pinned at zero so it stays out of the step-size
      // error norm -- see the note on aux_num.
      vars.set_rate(soil_number_of_depths + 4, 0.0);

  }

  // calculate K from K_sat based on theta
  double soil_K_from_soil_theta(double theta, size_t layer) const {
    //Eq. 5 Zeng and Decker (2009), ref Clapp and Hornberger (1978)
    // Clamp theta to [0, soil_moist_sat] (issue #485/#549): an intermediate
    // explicit-RK stage can probe theta < 0 (std::pow(negative, non-integer) is
    // NaN -> a non-positive layer drains nothing, K = 0) or, once the soil
    // feedback is perturbed, a theta far above saturation, whose large positive
    // exponent would otherwise make the inter-layer water_flux astronomically
    // large and cascade the blow-up across layers. Physically K saturates at
    // K_sat, so clamp there. Per-layer parameters (#558) fall back to the
    // scalar default when no layered vector is set.
    const double k_sat_layer = soil_parameter_value(K_sat_layers, K_sat, layer);
    const double soil_moist_sat_layer =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    const double n_psi_layer = soil_parameter_value(n_psi_layers, n_psi, layer);
    const double t = std::min(std::max(theta, 0.0), soil_moist_sat_layer);
    return k_sat_layer * std::pow(t / soil_moist_sat_layer, 2 * n_psi_layer + 3);
  }

  double soil_K_from_soil_theta(double theta) {
    return soil_K_from_soil_theta(theta, 0);
  }


  // convert soil moisture to soil water potential
  double psi_from_soil_moist(double soil_moist_, size_t layer) const {
    // Floor at the residual moisture: the retention curve (negative exponent)
    // diverges to +inf as theta->0, so an empty layer would otherwise yield a
    // non-finite potential. At/below theta_r the potential is large but finite
    // and the plant's root vulnerability curve has already shut uptake to ~0.
    const double t = std::max(soil_moist_, soil_moist_residual);
    const double a_psi_layer = soil_parameter_value(a_psi_layers, a_psi, layer);
    const double soil_moist_sat_layer =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    const double n_psi_layer = soil_parameter_value(n_psi_layers, n_psi, layer);
    const double psi =
      a_psi_layer * std::pow(t / soil_moist_sat_layer, -n_psi_layer) / 1e6; // Pa -> MPa
    // Cap at a large-but-finite ceiling (#549): near theta_r the retention curve
    // gives psi ~ 1e8 MPa, which pushes the leaf hydraulic solve non-finite and
    // corrupts the soil feedback. Root conductance is already ~0 far below this
    // (psi_crit ~ few MPa), so clamping to soil_psi_max_ leaves uptake at ~0
    // while keeping every downstream calculation finite.
    return std::min(psi, soil_psi_max_);
  }

  double psi_from_soil_moist(double soil_moist_) const {
    return psi_from_soil_moist(soil_moist_, 0);
  }

  // convert soil water potential to soil moisture
  double soil_moist_from_psi(double psi_soil_, size_t layer) const {
    const double a_psi_layer = soil_parameter_value(a_psi_layers, a_psi, layer);
    const double n_psi_layer = soil_parameter_value(n_psi_layers, n_psi, layer);
    const double soil_moist_sat_layer =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    // psi_soil_ is in MPa; a_psi is in Pa.
    return pow((psi_soil_ * 1e6 / a_psi_layer), (-1 / n_psi_layer)) * soil_moist_sat_layer;
  }

  double soil_moist_from_psi(double psi_soil_) const {
    return soil_moist_from_psi(psi_soil_, 0);
  }

  // Easy wrappers. Cn also use `extrinsic_drivers_evaluate("PPFD", time)
  //
  // Each is read once per individual per ODE derivs evaluation (in the leaf
  // physiology setup), always at the current `time`. Memoise per driver keyed
  // on `time` so the repeated unordered_map<string,...> lookups across
  // individuals at the same time collapse to one lookup per distinct time.
  // The memo is per-driver (not eager-refresh-all) so an unset driver is only
  // looked up if it is actually requested — preserving the existing throw-if-
  // absent behaviour. The returned value is bit-identical to a direct evaluate.
  double get_PPFD()       const { return cached_driver_("PPFD", ppfd_cache_, ppfd_cache_time_); }
  double get_atm_vpd()    const { return cached_driver_("atm_vpd", atm_vpd_cache_, atm_vpd_cache_time_); }
  double get_ca()         const { return cached_driver_("ca", ca_cache_, ca_cache_time_); }
  double get_leaf_temp()  const { return cached_driver_("leaf_temp", leaf_temp_cache_, leaf_temp_cache_time_); }
  double get_atm_o2_kpa() const { return cached_driver_("atm_o2_kpa", atm_o2_kpa_cache_, atm_o2_kpa_cache_time_); }
  double get_atm_kpa()    const { return cached_driver_("atm_kpa", atm_kpa_cache_, atm_kpa_cache_time_); }
  // Above-canopy wind speed U0 (m s^-1), Penman-Monteith aerodynamic resistance (#523)
  double get_wind_speed() const { return cached_driver_("wind_speed", wind_speed_cache_, wind_speed_cache_time_); }


  std::vector<double> get_soil_water_state() const { return {vars.states.begin(), vars.states.end() - aux_num}; }
  const std::vector<double>& get_soil_water_potential_state() const {
    bool cache_stale = !psi_soil_cache_valid_ ||
      psi_soil_cache_state_.size() != static_cast<size_t>(soil_number_of_depths);

    if (!cache_stale) {
      for (int i = 0; i < soil_number_of_depths; ++i) {
        if (psi_soil_cache_state_[i] != vars.state(i)) {
          cache_stale = true;
          break;
        }
      }
    }

    if (cache_stale) {
      psi_soil_cache_.resize(soil_number_of_depths);
      psi_soil_cache_state_.resize(soil_number_of_depths);
      for (int i = 0; i < soil_number_of_depths; ++i) {
        const double soil_moist = vars.state(i);
        psi_soil_cache_state_[i] = soil_moist;
        psi_soil_cache_[i] = psi_from_soil_moist(soil_moist, i);
      }
      psi_soil_cache_valid_ = true;
    }

    return psi_soil_cache_;
  }
  std::vector<double> get_soil_water_state_cumulative_flux() const { return {vars.states.end()-aux_num, vars.states.end()}; }
  std::vector<double> get_soil_depths() const { return z; }
  // double get_soil_depth(int layer) const { return z[layer]; }


  // TODO: I wonder if this needs a better name? See also environment.h
  Internals r_internals() const { return vars; }

  // The generic resource-pulse hook (#628): TF24's resources are its soil
  // layers, so a pulse is a depth of water added to one of them. Layer 0 is
  // rain; a deeper layer is irrigation, or a water table rising into the
  // profile. Either way the water has entered the column, so it counts as
  // input for the purposes of the water balance.
  std::vector<double> add_resource_pulse(size_t layer, double amount) {
    if (layer >= static_cast<size_t>(soil_number_of_depths)) {
      util::stop("Soil layer " + util::to_string(layer + 1) +
                 " does not exist: this environment has " +
                 util::to_string(soil_number_of_depths) + " layers");
    }
    return add_water_pulse_to_layer(layer, amount);
  }

  // Apply an instantaneous rainfall pulse of `depth` metres to the surface
  // layer. The model-specific name for the generic action above, kept because
  // rain onto the surface is what this is for nine times in ten.
  std::vector<double> add_water_pulse(double depth) {
    return add_resource_pulse(0, depth);
  }

  // R interface: layers are 1-based on that side, as everywhere else.
  std::vector<double> r_add_resource_pulse(util::index layer, double amount) {
    return add_resource_pulse(
      layer.check_bounds(static_cast<size_t>(soil_number_of_depths)), amount);
  }

  // The cap is the whole substance of this function. A pulse is applied
  // *between* solver legs, so unlike the continuous rates it has no error
  // estimate and no step rejection standing behind it: nothing but this line
  // stops it driving the layer past saturation, where the retention and
  // conductivity curves are meaningless. A layer can accept
  // (theta_sat - theta) * dz metres, and a realistic dryland event (~13 mm)
  // already exceeds that from a moderately wet start, so the excess is real
  // and frequent rather than a corner case. A hard min() is fine here
  // precisely because we are outside the integrator -- the "no kinks in the
  // rates" rule applies to continuous derivatives, not to a jump.
  //
  // Note this does *not* pass through the saturation-excess infiltration term
  // in compute_rates(): that term partitions a rainfall *rate* and is tuned
  // for continuous forcing, so applying it to an instantaneous depth would
  // double-count against the capacity cap. Whether a pulse should be filtered
  // that way as well is the first open question on this action.
  std::vector<double> add_water_pulse_to_layer(size_t layer, double depth) {
    if (!util::is_finite(depth) || depth < 0.0) {
      util::stop("Water pulse depth must be finite and non-negative");
    }
    const double theta = vars.state(layer);
    const double sat =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    const double capacity = std::max(0.0, (sat - theta) * dz[layer]);
    const double accepted = std::min(depth, capacity);
    const double excess = depth - accepted;

    vars.set_state(layer, theta + accepted / dz[layer]);
    // Direct state increments, not rates: the trailing slots are integrated
    // from their rates during a leg, and a pulse happens between legs.
    const size_t n = soil_number_of_depths;
    vars.set_state(n,     vars.state(n)     + depth);     // sum_rainfall
    vars.set_state(n + 1, vars.state(n + 1) + accepted);  // sum_infiltration
    vars.set_state(n + 4, vars.state(n + 4) + excess);    // sum_pulse_runoff

    psi_soil_cache_valid_ = false;
    return {accepted, excess};
  }

  // R interface
  void set_soil_water_state(std::vector<double> state) {
    if(state.size() != (vars.state_size- aux_num)) {
      throw std::invalid_argument("Input vector size does not match soil state size.");
    }
    for (size_t i = 0; i < (vars.state_size); i++) {
      if(i < soil_number_of_depths){
        vars.set_state(i, state[i]);
      } else {
        vars.set_state(i, 0);
      }
  }
    initial_states = vars.states;
    psi_soil_cache_valid_ = false;
}

  // Pre-compute resources available in the environment, as a function of height
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of light_availability spline
    // Note: extinction coefficient was already applied in strategy, so
    // f_compute_competition gives sum of projected leaf area (k L) across species. Just need to apply Beer's law, E = exp(- (k L))
    auto f_light_availability = [&](double height) -> double
    { return exp(-f_compute_competition(height)); };

    // Calculates the light_availability spline, by fitting to the function
    // `f_compute_competition` as a function of height
    light_availability.compute_environment(f_light_availability, height_max, rescale);
  }

  // Clear the light profile, and nothing else. This runs whenever the patch
  // empties (see StochasticPatch::compute_environment), not only between runs,
  // and the soil states are part of the ODE system: restoring them here would
  // discard what the solver had integrated, on every derivatives evaluation
  // while the patch held no individuals.
  virtual void clear_environment() {
    light_availability.clear();
  }

  // Return the soil moisture and the cumulative flux accumulators to the state
  // the run started from, so a second run on this patch starts where the first
  // did rather than silently continuing out of the first run's depleted soil.
  // Only Environment::clear() calls this.
  virtual void clear_state() {
    vars.states = initial_states;
    psi_soil_cache_valid_ = false;
  }

  virtual Rcpp::List r_get_state() const
  {
    
    // Surely an easier way?
    auto const &soil_depth_list = get_soil_depths();
    auto rcpp_soil_depth_vec = Rcpp::NumericVector(soil_depth_list.begin(), soil_depth_list.end());

    auto const &soil_moist_list = get_soil_water_state();
    auto rcpp_soil_moist_vec = Rcpp::NumericVector(soil_moist_list.begin(), soil_moist_list.end());

    auto const &soil_moist_cumulative_flux_list = get_soil_water_state_cumulative_flux();
    auto rcpp_soil_moist_vec_cumulative_flux = Rcpp::NumericVector(soil_moist_cumulative_flux_list.begin(), soil_moist_cumulative_flux_list.end());

    return Rcpp::List::create(
        // auto ret = get_state(environment.extrinsic_drivers, time);

        _["light_availability"] = light_availability.r_get_state(),
        _["soil_moist"] = rcpp_soil_moist_vec,
        _["soil_depth"] = rcpp_soil_depth_vec,
        _["soil_moist_cumulative_flux"] = rcpp_soil_moist_vec_cumulative_flux
    );
  }
  };
}

#endif
