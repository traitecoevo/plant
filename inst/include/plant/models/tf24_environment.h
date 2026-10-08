// Built from  inst/include/plant/models/ff16_environment.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_ENVIRONMENT_H_
#define PLANT_PLANT_TF24_ENVIRONMENT_H_

#include <plant/clamp_sites.h>
#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <array>
#include <limits>
#include <cmath>
#include <type_traits>

using namespace Rcpp;

namespace plant {

// Templated on the scalar S the light profile and the soil water balance both
// carry; double is production. The integrated soil state lives here rather than
// on the base because this is where its scalar is known: a cohort reads a layer
// potential derived from that state, so holding it at a value severs every
// route from a trait to the soil and back.
template <typename S = double>
class TF24_Environment : public Environment {
public:
  using value_type = S;

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
  
  TF24_Environment(int soil_number_of_depths = 5, 
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
    
    vars = Internals<S>(soil_number_of_depths + aux_num);
    initial_states = rebind_from(vars.states);

    z.resize(soil_number_of_depths);
    z_mid.resize(soil_number_of_depths);
    dz.resize(soil_number_of_depths);
    // positive downwards
    water_flux.resize(soil_number_of_depths);
    resize_resource_uptake();

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

    psi_soil_.resize(soil_number_of_depths);
    psi_soil_valid_ = false;

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
    psi_soil_valid_ = false;
  }

  int get_soil_number_of_depths() const {return soil_number_of_depths;}

  // The integrated soil state: one moisture per layer, then the four
  // cumulative-flux accumulators.
  size_t ode_size() const { return vars.state_size(); }

  template <typename It> It set_ode_state(It it) {
    for (size_t i = 0; i < vars.state_size(); i++) {
      vars.states[i] = *it++;
    }
    // The potentials are derived from the state just written, so what is held is
    // stale the moment this returns.
    psi_soil_valid_ = false;
    return it;
  }

  template <typename It> It ode_state(It it) const {
    for (size_t i = 0; i < vars.state_size(); i++) {
      util::write_iterator_scalar(it, vars.states[i]);
    }
    return it;
  }

  template <typename It> It ode_rates(It it) const {
    for (size_t i = 0; i < vars.state_size(); i++) {
      util::write_iterator_scalar(it, vars.rates[i]);
    }
    return it;
  }

  template <typename It> It ode_aux(It it) const {
    util::check_length(resource_uptake.size(), aux_size());
    for (size_t i = 0; i < aux_size(); i++) {
      util::write_iterator_scalar(it, resource_uptake[i]);
    }
    return it;
  }

  template <typename It> It set_ode_aux(It it) {
    util::check_length(resource_uptake.size(), aux_size());
    for (size_t i = 0; i < aux_size(); i++) {
      resource_uptake[i] = *it++;
    }
    return it;
  }

  // n_resources() is the count and this is the buffer it sizes, so an
  // environment that changes its resource count calls this and the two cannot
  // disagree.
  void resize_resource_uptake() {
    resource_uptake.assign(n_resources(), 0.0);
  }

  // The values of an active vector, for the R boundary and for the state a run
  // restarts from.
  static std::vector<double> rebind_from(const std::vector<S>& x) {
    std::vector<double> out(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      out[i] = odelia::util::to_passive(x[i]);
    }
    return out;
  }

  // Water is taken up per soil layer; the trailing aux_num state slots only
  // accumulate diagnostics and are never consumed.
  size_t n_resources() const override {
    return static_cast<size_t>(soil_number_of_depths);
  }

  // A cohort reads the light field and its own layers' water potentials. The
  // order is knot values, then knot slopes, then one potential per soil layer.
  // The count and the two fills below are declared together, so a reader that
  // gets the count from this class gets the fill from it too.
  size_t n_cohort_reads() const {
    return 2 * light_availability.knot_count() + 1 +
           static_cast<size_t>(soil_number_of_depths);
  }

  template <typename It> It cohort_reads(It it) const {
    const std::vector<S>& y = light_availability.knot_values();
    const std::vector<S>& m = light_availability.knot_slopes();
    util::check_length(y.size(), light_availability.knot_count());
    util::check_length(m.size(), light_availability.knot_count());
    for (size_t k = 0; k < y.size(); ++k) {
      util::write_iterator_scalar(it, y[k]);
    }
    for (size_t k = 0; k < m.size(); ++k) {
      util::write_iterator_scalar(it, m[k]);
    }
    // The canopy top, because the field is held against height / height_max and
    // a query divides by it: read the values and slopes without this and a
    // cohort's reads name every knot and not the grid they sit on.
    util::write_iterator_scalar(it, light_availability.height_max());
    const std::vector<S>& psi = get_soil_water_potential_state();
    for (int i = 0; i < soil_number_of_depths; ++i) {
      util::write_iterator_scalar(it, psi[i]);
    }
    return it;
  }

  template <typename It> It set_cohort_reads(It it) {
    const size_t n_knot = light_availability.knot_count();
    std::vector<S> y(n_knot), m(n_knot);
    for (size_t k = 0; k < n_knot; ++k) { y[k] = *it++; }
    for (size_t k = 0; k < n_knot; ++k) { m[k] = *it++; }
    light_availability.set_knot_data(y, m);
    light_availability.set_height_max(*it++);
    // ⚠️ MARKED VALID SO THE INJECTED POTENTIALS SURVIVE. They are not what the
    // state implies, and a read that derived them again would quietly replace
    // them -- which is the second thing the flag is for.
    psi_soil_.resize(soil_number_of_depths);
    for (int i = 0; i < soil_number_of_depths; ++i) {
      psi_soil_[i] = *it++;
    }
    psi_soil_valid_ = true;
    return it;
  }

  std::vector<double> get_soil_mid_depths() const { return z_mid; }

  // Another environment's values, written into this one.
  //
  // The light spline is not among them, and that is the same decision a rebind
  // makes: it is derived from the state, and whoever sets a state rebuilds it.
  // Its grid is 65 fixed fractions of the canopy top, so its count cannot drift
  // with what it holds, and its positions stay double while its values and slopes
  // carry the derivative.
  template <class S1>
  void assign_from(const TF24_Environment<S1>& src) {
    static_cast<Environment&>(*this) = static_cast<const Environment&>(src);
    vars = Internals<S>(src.vars.state_size(), src.vars.aux_size(),
                        src.vars.resource_size());
    for (size_t i = 0; i < src.vars.state_size(); ++i) {
      vars.states[i] = S(odelia::util::to_passive(src.vars.states[i]));
    }
    water_flux.assign(src.water_flux.size(), S(0.0));
    resource_uptake.assign(src.resource_uptake.size(), 0.0);
    z = src.z;
    z_mid = src.z_mid;
    dz = src.dz;
    initial_states = src.initial_states;
    soil_number_of_depths = src.soil_number_of_depths;
    delta_z = src.delta_z;
    depth = src.depth;
    soil_moist_sat = src.soil_moist_sat;
    K_sat = src.K_sat;
    a_psi = src.a_psi;
    n_psi = src.n_psi;
    soil_moist_sat_layers = src.soil_moist_sat_layers;
    K_sat_layers = src.K_sat_layers;
    a_psi_layers = src.a_psi_layers;
    n_psi_layers = src.n_psi_layers;
    use_layered_soil_parameters = src.use_layered_soil_parameters;
    a_infil = src.a_infil;
    b_infil = src.b_infil;
    soil_moist_residual = src.soil_moist_residual;
    soil_psi_max_ = src.soil_psi_max_;
    // Shared, not copied, for the reason the strategy's is: this copy is the one
    // the sweep runs on and it is discarded afterwards.
    clamps.differentiated = src.clamps.differentiated;
    // Derived from the state, so it is emptied rather than carried and the flag
    // says so.
    psi_soil_.assign(soil_number_of_depths, 0.0);
    psi_soil_valid_ = false;
    // The base assignment above carries `time`, so a memo left as it was would
    // read fresh at the new time holding the old environment's values.
    drivers_ = driver_memo{};
  }

  template <class U>
  TF24_Environment<U> rebind_from() const {
    TF24_Environment<U> out;
    out.assign_from(*this);
    return out;
  }

  // The state the solver integrates, and the uptake compute_rates received.
  Internals<S> vars;
  // Double: this is the uptake read back out through the R-facing aux interface,
  // never an input to a rate, so holding it at the working scalar was a slot per
  // layer per stage for a reading.
  std::vector<double> resource_uptake;

  // TODO: should we use auxilliary in internals
  std::vector<S> water_flux;
  std::vector<double> z;
  std::vector<double> z_mid;
  std::vector<double> dz;
  // The potentials the moisture state implies, derived where they are asked for
  // and kept until the state moves. Every writer of that state clears the flag
  // below, which is the whole of what keeps this true -- a second copy of the
  // state was kept beside these and compared per read, and it could only ever
  // catch a writer that had forgotten.
  //
  // ⚠️ A value read here cannot be compared against a value to decide staleness:
  // on an active scalar the derivative can move behind an unchanged value, and
  // the potentials would then carry the previous evaluation's tape slots.
  mutable std::vector<S> psi_soil_;
  mutable bool psi_soil_valid_ = false;

  static constexpr double NAN_TIME_ = std::numeric_limits<double>::quiet_NaN();

  // The time-varying drivers the leaf physiology reads. A driver is added here
  // and in one getter, and its name is read from this table rather than spelled
  // at the call site, so a memo cannot be keyed on one name and read under
  // another.
  enum driver {
    DRIVER_PPFD = 0,
    DRIVER_ATM_VPD,
    DRIVER_CA,
    DRIVER_LEAF_TEMP,
    DRIVER_ATM_O2_KPA,
    DRIVER_ATM_KPA,
    DRIVER_WIND_SPEED,
    DRIVER_COUNT
  };
  static constexpr std::array<const char*, DRIVER_COUNT> driver_names_ = {
      "PPFD", "atm_vpd", "ca", "leaf_temp", "atm_o2_kpa", "atm_kpa",
      "wind_speed"};

  // What was read, and the one time every entry was read at. Each driver is a
  // function of that time and of nothing else, so the memo holds one time and
  // says per driver whether it has been asked for yet.
  struct driver_memo {
    double time = NAN_TIME_;
    std::array<bool, DRIVER_COUNT> read{};
    std::array<double, DRIVER_COUNT> value{};
  };
  mutable driver_memo drivers_;

  // ⚠️ ONE DRIVER AT A TIME, AND NOT ALL SEVEN WHERE THE TIME IS SET. `evaluate`
  // raises for a driver that was never set and for a time outside a variable
  // driver's control points, so refreshing seven to serve one would raise on an
  // environment that only ever reads one of them.
  double driver_at(int d) const {
    if (drivers_.time != time) {
      drivers_.time = time;
      drivers_.read.fill(false);
    }
    if (!drivers_.read[d]) {
      drivers_.value[d] = extrinsic_drivers.evaluate(driver_names_[d], time);
      drivers_.read[d] = true;
    }
    return drivers_.value[d];
  }

  // The soil state a run begins from: whatever set_soil_water_state last set,
  // which the R interface lets a caller choose. Restored by clear_state().
  std::vector<double> initial_states;

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline<S> light_availability;

  // Every active value this environment holds: the integrated state and its
  // rates, the per-layer flux, the potentials derived from the state, and the
  // light field's knot values and slopes.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, vars, water_flux, psi_soil_,
                              light_availability);
  }

  // Light interface
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

  // The clamp sites the water balance reaches. The site list is shared with the
  // strategy's, so a tally reads the same however the site is reached; the
  // differentiated half is shared through rebind_from for the same reason.
  // Mutable because it is a tally and not state: the reads that clamp are const,
  // and a clamp changes what they return to nobody.
  mutable clamp_counter clamps;
  void note_clamp(int site) const {
    if constexpr (std::is_same_v<S, double>) {
      ++clamps.forward[site];
    } else {
      ++(*clamps.differentiated)[site];
    }
  }

  // Ability to prescribe a fixed value
  // TODO: add setting to set other variables like water
  void set_fixed_environment(S value, S height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(S value) {
    S height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  S get_environment_at_height(S height) const {
    return light_availability.get_value_at_height(height);
  }

  // Highest height covered by the light spline; hoist out of hot per-point
  // loops and feed back into the capped get_environment_at_height() overload.
  double max_environment_height() const {
    return light_availability.max_height();
  }

  S get_environment_at_height(S height, double cap) const {
    return light_availability.get_value_at_height(height, cap);
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }
  virtual std::vector<double> get_interpolators_state() const {
    return light_availability.interpolators_state();
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
  void compute_rates(std::vector<S> const &resource_depletion)
  {
    using std::max;
    using std::pow;

    S water_input;
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
    const double rainfall_raw = extrinsic_drivers.evaluate("rainfall", time);
    if (rainfall_raw < 0.0) {
      note_clamp(CLAMP_RAINFALL);
    }
    double rainfall = std::max(0.0, rainfall_raw);
    const double soil_moist_sat_0 =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, 0);
    const S excess = S(1) - a_infil * pow(vars.state(0) / soil_moist_sat_0, b_infil);
    // Counted: where the top layer is wet enough that this goes negative, the
    // infiltrated fraction stops reading that layer's state at all.
    if (excess < S(0.0)) {
      note_clamp(CLAMP_INFILTRATION);
    }
    S infiltration = rainfall * max(S(0.0), excess);
    S total_resource_depletion = 0.0;


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
      const S theta = vars.state(i);
      S rate = (water_input - water_flux[i] - resource_depletion[i]) / dz[i];
      // Positivity guard (issue #485), hardened for #549: at/below the residual
      // moisture a layer must not be dried further. The original `rate < 0.0`
      // test let a *non-finite* rate through, because `NaN < 0.0` is false in
      // IEEE 754 -- a NaN soil_consumption_ (from the retention curve's enormous
      // psi_soil near theta_r) then wrote straight into the soil state. `!(rate
      // > 0.0)` is true for NaN and for rate <= 0, so only genuine rewetting is
      // allowed and NaN/negative rates are clamped to 0.
      if (theta <= soil_moist_residual && !(rate > 0.0)) {
        // Counted: a layer held here reads no uptake and no flux, so every route
        // from a trait to this layer's state is cut for as long as it holds.
        note_clamp(CLAMP_SOIL_POSITIVITY);
        rate = 0.0;
      }
      vars.set_rate(i, rate);
      resource_uptake[i] = odelia::util::to_passive(resource_depletion[i]);
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
  S soil_K_from_soil_theta(S theta, size_t layer) const {
    using std::max;
    using std::min;
    using std::pow;
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
    if (theta < S(0.0) || theta > S(soil_moist_sat_layer)) {
      note_clamp(CLAMP_SOIL_CONDUCTIVITY);
    }
    const S t = min(max(theta, S(0.0)), S(soil_moist_sat_layer));
    return k_sat_layer * pow(t / soil_moist_sat_layer, 2 * n_psi_layer + 3);
  }

  S soil_K_from_soil_theta(S theta) {
    return soil_K_from_soil_theta(theta, 0);
  }


  // convert soil moisture to soil water potential
  S psi_from_soil_moist(S soil_moist_, size_t layer) const {
    using std::max;
    using std::min;
    using std::pow;
    // Floor at the residual moisture: the retention curve (negative exponent)
    // diverges to +inf as theta->0, so an empty layer would otherwise yield a
    // non-finite potential. At/below theta_r the potential is large but finite
    // and the plant's root vulnerability curve has already shut uptake to ~0.
    if (soil_moist_ < S(soil_moist_residual)) {
      note_clamp(CLAMP_SOIL_MOISTURE_FLOOR);
    }
    const S t = max(soil_moist_, S(soil_moist_residual));
    const double a_psi_layer = soil_parameter_value(a_psi_layers, a_psi, layer);
    const double soil_moist_sat_layer =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    const double n_psi_layer = soil_parameter_value(n_psi_layers, n_psi, layer);
    const S psi =
      a_psi_layer * pow(t / soil_moist_sat_layer, -n_psi_layer) / 1e6; // Pa -> MPa
    // Cap at a large-but-finite ceiling (#549): near theta_r the retention curve
    // gives psi ~ 1e8 MPa, which pushes the leaf hydraulic solve non-finite and
    // corrupts the soil feedback. Root conductance is already ~0 far below this
    // (psi_crit ~ few MPa), so clamping to soil_psi_max_ leaves uptake at ~0
    // while keeping every downstream calculation finite.
    if (psi > S(soil_psi_max_)) {
      note_clamp(CLAMP_SOIL_POTENTIAL_CEILING);
    }
    return min(psi, S(soil_psi_max_));
  }

  S psi_from_soil_moist(S soil_moist_) const {
    return psi_from_soil_moist(soil_moist_, 0);
  }

  // convert soil water potential to soil moisture
  double soil_moist_from_psi(double psi, size_t layer) const {
    const double a_psi_layer = soil_parameter_value(a_psi_layers, a_psi, layer);
    const double n_psi_layer = soil_parameter_value(n_psi_layers, n_psi, layer);
    const double soil_moist_sat_layer =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    // psi is in MPa; a_psi is in Pa.
    return pow((psi * 1e6 / a_psi_layer), (-1 / n_psi_layer)) * soil_moist_sat_layer;
  }

  double soil_moist_from_psi(double psi) const {
    return soil_moist_from_psi(psi, 0);
  }

  // Each is read once per individual per ODE derivs evaluation, always at the
  // current `time`, so the memo collapses the repeated string-keyed lookups
  // across individuals into one lookup per driver per distinct time. Every one
  // returns what a direct `extrinsic_drivers_evaluate(name, time)` returns.
  double get_PPFD()       const { return driver_at(DRIVER_PPFD); }
  double get_atm_vpd()    const { return driver_at(DRIVER_ATM_VPD); }
  double get_ca()         const { return driver_at(DRIVER_CA); }
  double get_leaf_temp()  const { return driver_at(DRIVER_LEAF_TEMP); }
  double get_atm_o2_kpa() const { return driver_at(DRIVER_ATM_O2_KPA); }
  double get_atm_kpa()    const { return driver_at(DRIVER_ATM_KPA); }
  // Above-canopy wind speed U0 (m s^-1), Penman-Monteith aerodynamic resistance (#523)
  double get_wind_speed() const { return driver_at(DRIVER_WIND_SPEED); }


  std::vector<double> get_soil_water_state() const {
    std::vector<double> out = rebind_from(vars.states);
    out.resize(out.size() - aux_num);
    return out;
  }
  const std::vector<S>& get_soil_water_potential_state() const {
    if (!psi_soil_valid_ ||
        psi_soil_.size() != static_cast<size_t>(soil_number_of_depths)) {
      psi_soil_.resize(soil_number_of_depths);
      for (int i = 0; i < soil_number_of_depths; ++i) {
        psi_soil_[i] = psi_from_soil_moist(vars.state(i), i);
      }
      psi_soil_valid_ = true;
    }
    return psi_soil_;
  }
  std::vector<double> get_soil_water_state_cumulative_flux() const {
    const std::vector<double> all = rebind_from(vars.states);
    return {all.end() - aux_num, all.end()};
  }
  std::vector<double> get_soil_depths() const { return z; }
  // double get_soil_depth(int layer) const { return z[layer]; }


  // TODO: I wonder if this needs a better name? See also environment.h
  Internals<double> r_internals() const {
    Internals<double> out(vars.state_size(), vars.aux_size(), vars.resource_size());
    out.states = rebind_from(vars.states);
    out.rates = rebind_from(vars.rates);
    return out;
  }

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
    // ⚠️ THE STATE IS AT THIS ENVIRONMENT'S SCALAR, NOT double. A pulse is an
    // instantaneous state increment applied between legs, so it is arithmetic on
    // whatever the state is -- and the two doubles crossing the boundary are the
    // ARGUMENT and the report, which is why they are read passive here rather
    // than the state being narrowed.
    const S theta = vars.state(layer);
    const S sat = S(
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer));
    const S capacity_raw = (sat - theta) * dz[layer];
    const double capacity =
        std::max(0.0, odelia::util::to_passive(capacity_raw));
    const double accepted = std::min(depth, capacity);
    const double excess = depth - accepted;

    vars.set_state(layer, theta + S(accepted / dz[layer]));
    // Direct state increments, not rates: the trailing slots are integrated
    // from their rates during a leg, and a pulse happens between legs.
    const size_t n = soil_number_of_depths;
    vars.set_state(n,     vars.state(n)     + S(depth));     // sum_rainfall
    vars.set_state(n + 1, vars.state(n + 1) + S(accepted));  // sum_infiltration
    vars.set_state(n + 4, vars.state(n + 4) + S(excess));    // sum_pulse_runoff

    psi_soil_valid_ = false;
    return {accepted, excess};
  }

  // R interface
  void set_soil_water_state(std::vector<double> state) {
    if(state.size() != (vars.state_size()- aux_num)) {
      throw std::invalid_argument("Input vector size does not match soil state size.");
    }
    for (size_t i = 0; i < (vars.state_size()); i++) {
      if(i < soil_number_of_depths){
        vars.set_state(i, state[i]);
      } else {
        vars.set_state(i, 0);
      }
  }
    initial_states = rebind_from(vars.states);
    psi_soil_valid_ = false;
}

  // The light a height is left with, from the competition profile above it.
  template <typename Function>
  void compute_environment(Function f_compute_competition_and_slope,
                           S height_max) {
    build_extinction_field(light_availability, f_compute_competition_and_slope,
                           height_max);
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
    util::check_length(initial_states.size(), vars.state_size());
    for (size_t i = 0; i < vars.state_size(); ++i) {
      vars.states[i] = S(initial_states[i]);
    }
    psi_soil_valid_ = false;
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
