// -*-c++-*-
#ifndef PLANT_PLANT_TF24_ENVIRONMENT_H_
#define PLANT_PLANT_TF24_ENVIRONMENT_H_

#include <plant/gradients/clamp_sites.h>
#include <plant/models/environment.h>
#include <plant/models/resource_spline.h>
#include <plant/models/leaf_model.h>
#include <algorithm>
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
class TF24_Environment final : public Environment {
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

  // constructor for R interface - default settings can be modified.
  // Layer geometry, soil_number_of_depths included, can be changed after
  // construction, through set_soil_number_of_depths or set_soil_layer_widths.
  TF24_Environment(int soil_number_of_depths = 5,
                   double soil_moist_sat = 0.428, // saturated soil moisture content (m3 water m^-3 soil)
                   double K_sat = 163.0411, //saturated hydraulic conductivity of soil
                   double a_psi = 1.78e-3, // retention-curve scale (MPa)
                   double n_psi = 6.57, // retention-curve exponent (unitless)
                   double a_infil = 1, // infiltration switch (0-1), 0 no runoff, 1 runoff
                   double b_infil = 8, // unitless, determines infiltration rate
                   double depth = 1.5)  // total depth of soil (m)
      : soil_moist_sat(soil_moist_sat),
      a_psi(a_psi),
      n_psi(n_psi),
      K_sat(K_sat),
      a_infil(a_infil),
      b_infil(b_infil),
      depth(depth)
  {
    time = 0.0;

    // The light driver is incoming SHORTWAVE at the surface, on a horizontal plane
    // above the canopy (W m^-2) -- the quantity a flux tower reports as SWdown, and
    // NOT a top-of-atmosphere value. 900 is a clear-sky midday figure. The
    // strategy splits it into its PAR and near-infrared bands with the f_PAR
    // parameter and carries each down to the leaf; see
    // TF24_Model::compute_leaf_physiology.
    extrinsic_drivers_set_constant("SWdown",900);
    extrinsic_drivers_set_constant("rainfall",1);
    extrinsic_drivers_set_constant("atm_vpd",1);
    extrinsic_drivers_set_constant("ca",40);
    extrinsic_drivers_set_constant("leaf_temp",25);
    extrinsic_drivers_set_constant("atm_o2_kpa",21);
    // 101.3 kPa, standard sea-level pressure, and the value every leaf-level
    // test uses. The leaf derives its ppm -> Pa conversion from atm_kpa, so the
    // model is self-consistent at any pressure, but the published numbers are at
    // this one: 100.5 kPa shifts TF24 output by 2.4%. It is not a site elevation;
    // set it per-site if you mean altitude.
    extrinsic_drivers_set_constant("atm_kpa",101.3);
    extrinsic_drivers_set_constant("wind_speed",2.0); // U0, m s^-1 (Penman-Monteith)

    set_soil_number_of_depths(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, soil_moist_sat*0.5));
  };
  
  // Number of cumulative auxilliary variables to track in soil moisture model
  // Trailing state slots that accumulate diagnostic fluxes rather than
  // participating in the water balance: sum_rainfall, sum_infiltration,
  // sum_drainage, sum_resource_depletion, sum_pulse_runoff.
  //
  // The fifth is fed only by rainfall-pulse events. That is deliberate,
  // not an oversight: the ODE controller takes its error norm over *every*
  // state component, so an accumulator with a non-zero rate would join the
  // step-size decision and could move existing TF24 trajectories. Held at rate
  // zero it contributes exactly zero to the norm, so a run with no pulses is
  // bit-identical to one from before this slot existed. Accumulating the
  // continuous saturation-excess runoff here too would be a deliberate change
  // to what TF24 reports, and needs its own decision.
  static constexpr size_t aux_num = 5;
  
  // ------------------------------------------------------------------
  // SOIL LAYER GEOMETRY
  // ------------------------------------------------------------------
  // Four descriptions of the same discretisation, and each is read by a
  // different part of the model, so all four are kept rather than derived on
  // demand:
  //
  //   dz[i]    the WIDTH of layer i (m). The water balance divides by it, and the
  //            pulse capacity cap multiplies by it. This is the geometry a caller
  //            states; the other two follow from it.
  //   z[i]     the cumulative depth to the BOTTOM of layer i, so z.back() is the
  //            total column depth. Read by the root-mass distribution, which
  //            differences Q(z) across a layer -- correct only because these are
  //            boundaries and not sample depths.
  //   z_mid[i] the MIDPOINT of layer i. Read by the leaf's gravitational head.
  //   z_thickness[i] the width of layer i as z's differences, through
  //            phylloptim's layer_thicknesses. Read by the root network. It
  //            equals dz only to rounding, so neither stands in for the other.
  //
  // ⚠️ WIDTH IS THE INPUT, NOT DEPTH. The two setters below differ only in
  // which they take, and a caller wanting unequal layers must go through
  // set_soil_layer_widths: there is deliberately no way to state a boundary
  // directly, because "depth 0.3" is ambiguous between a top, a midpoint and a
  // bottom while "width 0.3" is not.

  // Equal layers: n of them, spanning the current `depth`.
  void set_soil_number_of_depths(int n) {
    soil_number_of_depths = n;

    z.resize(soil_number_of_depths);
    dz.resize(soil_number_of_depths);

    // TODO(next rebaseline): build `z` as the running sum of `dz`, as
    // set_soil_layer_widths does. `(i+1)*delta_z` and that sum disagree in the
    // last bit for 103 of 180 (depth, n) pairs, and `z` feeds Q() and the
    // gravitational head, so the change moves every TF24 number.
    const double delta_z = depth / soil_number_of_depths;
    for (int i = 0; i < soil_number_of_depths; i++)
    {
      z[i] = (i + 1) * delta_z;
      dz[i] = delta_z;
    }

    finalise_soil_geometry_();
  }

  // Unequal layers: the width of each, top down. `depth` becomes their total.
  //
  // This is what soil evaporation needs -- a very thin layer at the top of the
  // profile, which cannot be expressed as a count.
  //
  // ⚠️ A THIN SURFACE LAYER MAKES THE COLUMN STIFF, and nothing here stops you.
  // The water balance's diagonal is K'(theta)/dz, so it scales as 1/dz: measured
  // 1.8e4 yr^-1 at the 0.3 m default, 2.6e5 at 2 cm, where the characteristic time
  // is about two minutes. It only bites when the layer is wet, i.e. just after
  // rain, which is exactly when an evaporation layer matters. Below dz ~ 2.6 mm
  // the stable step falls under Control()'s ode_step_size_min and the controller
  // cannot satisfy stability at all. See .claude/hazards.md, hazard 11.
  //
  // ⚠️ AND IT RECALIBRATES RUNOFF. The infiltration term keys off layer 0 alone
  // through (theta_0/theta_sat)^b_infil with b_infil = 8, so a thin layer 0 fills
  // far faster and holds the switch shut for longer. a_infil/b_infil were
  // calibrated against a 30 cm surface layer. Not corrected here, because doing so
  // is a change to what TF24 reports and needs its own decision.
  void set_soil_layer_widths(std::vector<double> widths) {
    if (widths.empty()) {
      util::stop("set_soil_layer_widths: need at least one layer");
    }
    for (size_t i = 0; i < widths.size(); ++i) {
      if (!util::is_finite(widths[i]) || widths[i] <= 0.0) {
        util::stop("set_soil_layer_widths: every width must be finite and "
                   "positive; layer " + util::to_string(i + 1) + " is " +
                   util::format_double(widths[i]));
      }
    }
    // Layers must not get THINNER with depth. Equal is fine, so a uniform
    // profile is legal and `rep(w, n)` still works.
    //
    // WHY. Refining towards the surface is the standard: a 1 cm-scale cell is
    // needed near the surface and not deeper down (Downer & Ogden 2004), and
    // that is what the profile here is for. Refining at DEPTH has no counterpart
    // in this model -- a contrast in soil properties is per-layer *parameters*
    // (set_soil_parameters), not a thin layer, and a water table is a
    // bottom boundary condition TF24 does not have -- while it does reach a
    // failure this model cannot absorb. A thin layer under a thick one drives
    // the size-density equation to overflow: sampled over the basal
    // thickness, it fails at 5 cm, 3 cm, 2.5 cm, 2 cm and 1.5 cm, where all 18
    // profiles thickening with depth ran clean, including a 5 mm surface layer.
    //
    // ⚠️ THIS IS NOT A SAFETY GUARANTEE, and do not read it as one. The failures
    // it blocks are RAGGED in the geometry -- 5 cm fails, 4 cm passes, 3 to
    // 1.5 cm fail, 1 cm passes again -- so thinning-with-depth is a property the
    // observed failures share rather than a characterisation of them, and 46 of
    // 48 sampled thinning profiles ran fine. It is a foot-gun removed, not a
    // proof of stability; the underlying fragility is in the size-density
    // equations.
    //
    // ⚠️ NOR IS A THICKENING PROFILE AUTOMATICALLY WELL RESOLVED. The error in
    // the usual Richards discretisation vanishes on a uniform grid and "worsens
    // with increasing layer thickness differences" whichever way they run
    // (Mackay et al. 2022), and a fast increase with depth is equivalent to a
    // coarser column (Regenass et al. 2021). So grading steeply is its own cost;
    // see soil_widths_graded().
    for (size_t i = 1; i < widths.size(); ++i) {
      if (widths[i] < widths[i - 1]) {
        util::stop(
            "set_soil_layer_widths: layers must not get thinner with depth; "
            "layer " + util::to_string(i + 1) + " (" +
            util::format_double(widths[i]) + " m) is thinner than layer " +
            util::to_string(i) + " (" + util::format_double(widths[i - 1]) +
            " m). Equal widths are allowed. Refine towards the SURFACE, which "
            "is what a thin layer is for here; a thin layer beneath a thicker "
            "one can drive the size-density equations to overflow, and "
            "a contrast in soil properties belongs in set_soil_parameters() "
            "rather than in the geometry.");
      }
    }

    soil_number_of_depths = static_cast<int>(widths.size());
    dz = widths;

    z.resize(soil_number_of_depths);
    double cumulative = 0.0;
    for (int i = 0; i < soil_number_of_depths; i++) {
      cumulative += dz[i];
      z[i] = cumulative;
    }
    // Taken off the running sum rather than summed again, so `depth` and
    // `z.back()` cannot disagree by rounding. Everything downstream reads one or
    // the other and they must be the same number.
    depth = z.back();

    finalise_soil_geometry_();
  }

  // Everything that follows from z and dz, whichever setter produced them.
  //
  // ⚠️ THIS RESETS THE WATER STATE AND DROPS ANY LAYERED PARAMETERS. Both are
  // per-layer, so neither survives a change of layer count, and the geometry
  // setters are not the place to guess a re-mapping. Set the geometry first, then
  // the parameters, then the state -- or use plant's set_tf24_soil() from R, which
  // exists so nobody has to remember that order.
  void finalise_soil_geometry_() {
    vars = Internals<S>(soil_number_of_depths + aux_num);
    initial_states = rebind(vars.states);

    z_mid.resize(soil_number_of_depths);
    for (int i = 0; i < soil_number_of_depths; i++) {
      if (i == 0) {
        z_mid[i] = z[i] / 2.0;
      } else {
        z_mid[i] = (z[i - 1] + z[i]) / 2.0;
      }
    }
    z_thickness = phylloptim::layer_thicknesses(z);

    // positive downwards
    water_flux.resize(soil_number_of_depths);
    resize_resource_uptake();

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
    // Rebuild the geometry only if the layer count actually changed. Called after
    // set_soil_layer_widths() at the same count, an unconditional rebuild would
    // silently throw the widths away and revert to equal layers -- which is a
    // failure mode with no symptom, since both are valid geometries.
    if (n != soil_number_of_depths) {
      set_soil_number_of_depths(n);
    }

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

  // The integrated soil state: one moisture per layer, then the `aux_num`
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
  static std::vector<double> rebind(const std::vector<S>& x) {
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
    // The canopy top, because the field is tabulated against height / height_max and
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

  const std::vector<double>& get_soil_mid_depths() const { return z_mid; }

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
    resource_uptake.assign(src.resource_uptake.size(), S(0.0));
    z = src.z;
    z_mid = src.z_mid;
    z_thickness = src.z_thickness;
    dz = src.dz;
    initial_states = src.initial_states;
    soil_number_of_depths = src.soil_number_of_depths;
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
    f_PAR = src.f_PAR;
    soil_moist_residual = src.soil_moist_residual;
    soil_psi_max_ = src.soil_psi_max_;
    water_mixing_rate = src.water_mixing_rate;
    private_water_mean_ = src.private_water_mean_;
    private_water_newborn_ = src.private_water_newborn_;
    private_water_cover_ = src.private_water_cover_;
    // Shared, not copied, for the reason the strategy's is: this copy is the one
    // the sweep runs on and it is discarded afterwards.
    clamps.differentiated = src.clamps.differentiated;
    // And the light field's, which is rebuilt rather than carried.
    light_availability.clamps.differentiated =
        src.light_availability.clamps.differentiated;
    // Derived from the state, so it is emptied rather than carried and the flag
    // says so.
    psi_soil_.assign(soil_number_of_depths, 0.0);
    psi_soil_valid_ = false;
    // The base assignment above carries `time`, so a memo left as it was would
    // read fresh at the new time holding the old environment's values.
    drivers_ = driver_memo{};
  }

  template <class U>
  TF24_Environment<U> rebind() const {
    TF24_Environment<U> out;
    out.assign_from(*this);
    return out;
  }

  // The state the solver integrates, and the uptake compute_rates received.
  Internals<S> vars;

  // * Private root-zone water (TF24w)
  //
  // Lateral exchange between a plant's root-zone column and its surroundings,
  // as a relaxation rate toward the area-weighted patch mean [yr^-1]. A soil
  // property rather than a plant's, which is why it lives beside K_sat: one
  // number per patch, 0 is fully private, large is pooled (the TF24 limit).
  // Read by TF24w and by Patch's exchange pass; TF24 reads nothing from it.
  double water_mixing_rate = 0.0;
  // The aggregates the exchange pass computes from every column before the
  // species are rated, and that every column then reads: the area-weighted
  // mean water content per layer, the crown cover (sum of crown area per m2 of
  // patch) and the bare fraction max(0, 1 - cover). With no plants the mean is
  // this column, so a newborn seeded from it starts where the soil is.
  std::vector<double> private_water_mean_;
  // What a newborn's column starts at: the content of the ground it takes over,
  // which is the bare column while ground is open and the private columns' mean
  // once crowns cover it. Seeding at that content is what keeps an introduction
  // from making or losing water.
  std::vector<double> private_water_newborn_;
  double private_water_cover_ = 0.0;
  void set_private_water_aggregates(const std::vector<double>& mean,
                                    const std::vector<double>& newborn,
                                    double cover) {
    private_water_mean_ = mean;
    private_water_newborn_ = newborn;
    private_water_cover_ = cover;
  }
  double private_water_newborn(size_t layer) const {
    if (private_water_newborn_.size() != static_cast<size_t>(soil_number_of_depths)) {
      return odelia::util::to_passive(vars.state(layer));
    }
    return private_water_newborn_[layer];
  }
  double private_water_mean(size_t layer) const {
    if (private_water_mean_.size() != static_cast<size_t>(soil_number_of_depths)) {
      return odelia::util::to_passive(vars.state(layer));
    }
    return private_water_mean_[layer];
  }
  double private_water_cover() const { return private_water_cover_; }
  double private_water_bare() const { return std::max(0.0, 1.0 - private_water_cover_); }
  // The exchange pass adds the mixing and land-transfer terms to this column's
  // rates after compute_rates has set its own, and overwrites the diagnostic
  // accumulators with area-weighted totals over every column.
  void add_to_rate(size_t layer, const S& v) {
    vars.set_rate(layer, vars.rate(layer) + v);
  }
  // k: 0 rainfall, 1 infiltration, 2 drainage, 3 resource depletion, 4 pulse runoff.
  void set_accumulator_rate(size_t k, const S& v) {
    vars.set_rate(soil_number_of_depths + k, v);
  }
  const S& accumulator_rate(size_t k) const {
    return vars.rate(soil_number_of_depths + k);
  }
  // A pulse applied to every column: this column's add_water_pulse_to_layer has
  // already booked its own accepted depth and excess as if it were the whole
  // patch, so the pass corrects the two accumulators by the area-weighted
  // totals minus what was booked.
  void adjust_pulse_accumulators(double d_infiltration, double d_runoff) {
    const size_t n = soil_number_of_depths;
    vars.set_state(n + 1, vars.state(n + 1) + S(d_infiltration));
    vars.set_state(n + 4, vars.state(n + 4) + S(d_runoff));
  }
  // Why a state the solver proposes is outside the soil's domain, read from
  // layer 0's entry at `first`: empty when every layer holds no more than its
  // saturation, and otherwise the first layer that holds more. The accumulators
  // after the layers are not bounded. Patch::ode_state_valid reads it, so a
  // step that ends above saturation is rejected and retried smaller rather
  // than committed: the moisture is a share of the pore space, and above
  // saturation the retention curve gives a suction smaller than a saturated
  // soil's. Asked as "is it inside", so a NaN fails.
  //
  // The balance itself keeps every layer at or below saturation (see
  // accepted_inflow), so a refusal is an explicit step's overshoot, which a
  // smaller step removes. The one inflow the balance does not limit is water
  // the roots release into a layer (a negative uptake): released into a full
  // layer it overfills it, no smaller step escapes that, and the run stops,
  // naming the layer.
  template <class It>
  std::string saturation_refusal(It first) const {
    // The bound allows a few ulps above saturation. A capped rain pulse fills
    // a layer to theta + (sat - theta) * dz / dz, which rounds to within an
    // ulp or two of sat, either side, depending on how the compiler contracts
    // it; and events are applied between steps, where nothing checks them. A
    // full layer with no net flux would otherwise hold a rounding-level excess
    // no step can remove, and every later step would be refused.
    constexpr double rounding = 16 * std::numeric_limits<double>::epsilon();
    for (int i = 0; i < soil_number_of_depths; ++i, ++first) {
      const double sat =
        soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, i);
      if (!(*first <= sat * (1.0 + rounding))) {
        return "soil layer " + util::to_string(i + 1) + " holds moisture " +
               util::format_double(*first) + ", above its saturation " +
               util::format_double(sat);
      }
    }
    return std::string();
  }
  // What a layer at `theta` can still accept, in metres of water: the cap
  // add_water_pulse_to_layer applies, factored out so a private column applies
  // the same one.
  double pulse_capacity(const S& theta, size_t layer) const {
    return std::max(0.0, odelia::util::to_passive(pulse_room(theta, layer)));
  }
  // The room below saturation at the state's scalar, unclamped: what a pulse
  // that fills the layer adds to it.
  S pulse_room(const S& theta, size_t layer) const {
    const S sat = S(
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer));
    return (sat - theta) * dz[layer];
  }
  // At the working scalar: a step that substeps the soil reads it as the
  // layers' inputs, and the sweep differentiates through them.
  std::vector<S> resource_uptake;

  std::vector<S> water_flux;
  std::vector<double> z;
  std::vector<double> z_mid;
  std::vector<double> z_thickness;
  std::vector<double> dz;
  // The potentials the moisture state implies, derived where they are asked for
  // and kept until the state moves. Every writer of that state clears the flag
  // below, which is the whole of what keeps this true -- a second copy of the
  // state kept beside these and compared per read could only ever catch a
  // writer that had forgotten.
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
    DRIVER_SWDOWN = 0,
    DRIVER_ATM_VPD,
    DRIVER_CA,
    DRIVER_LEAF_TEMP,
    DRIVER_ATM_O2_KPA,
    DRIVER_ATM_KPA,
    DRIVER_WIND_SPEED,
    DRIVER_COUNT
  };
  static constexpr std::array<const char*, DRIVER_COUNT> driver_names_ = {
      "SWdown", "atm_vpd", "ca", "leaf_temp", "atm_o2_kpa", "atm_kpa",
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
  // rates, the per-layer flux and uptake, the potentials derived from the state,
  // and the light field's knot values and slopes.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, vars, water_flux, resource_uptake, psi_soil_,
                              light_availability);
  }

  int soil_number_of_depths;

  // Total column depth, m: the sum of the layer widths, and equal to z.back() by
  // construction. READ-ONLY from R: only the geometry setters rebuild z/z_mid/z_thickness/dz,
  // so a raw write would leave them describing the old column. To change it, call set_soil_number_of_depths() again (equal layers spanning the
  // new depth) or set_soil_layer_widths() (which sets it from the widths).
  double depth;
  //saturated soil moisture
  double soil_moist_sat;
  //Saturated soil hydraulic conductivity
  double K_sat;
  // Retention curve psi = a_psi (theta / theta_sat)^-n_psi, a_psi in MPa.
  double a_psi;
  double n_psi;
  std::vector<double> soil_moist_sat_layers;
  std::vector<double> K_sat_layers;
  std::vector<double> a_psi_layers;
  std::vector<double> n_psi_layers;
  bool use_layered_soil_parameters = false;
  double a_infil;
  double b_infil;
  // PAR share of the incoming shortwave, on an ENERGY basis (unitless). A
  // property of the sky, about 0.45-0.5 and rising a little under cloud, so in
  // practice a constant: a fixed parameter rather than a driver.
  double f_PAR = 0.5;

  // Residual soil moisture floor (m3 m^-3). The water balance cannot dry a
  // layer below this (see the positivity guard in compute_rates), and
  // psi_from_soil_moist is evaluated at this floor so the retention curve
  // (which diverges as theta->0) stays finite. Well below any realistic
  // operating moisture, so it does not perturb non-drought runs.
  double soil_moist_residual = 1e-2;

  // Ceiling on soil water potential (MPa). Near the residual floor
  // the retention curve returns psi ~ 1e8 MPa, which drives the leaf hydraulic
  // solve non-finite and blows up the soil feedback. phylloptim bounds the root
  // vulnerability curves past their grid, so root
  // conductivity is exactly 0 above ~46 MPa at the default root_P50 = 3.4, a
  // point that scales with root_P50; at 50 MPa it is below 1e-22 for any
  // root_P50 up to 10. So the ceiling leaves uptake unchanged and tightening it
  // cannot fix a drought symptom. Dry runs do reach it, so it is set
  // near that point rather than at an oven-dry potential.
  double soil_psi_max_ = 50;

  // The clamp sites the water balance reaches. The site list is shared with the
  // strategy's, so a tally reads the same however the site is reached; the
  // differentiated half is shared through rebind for the same reason.
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
  //   * water_in_0   = infiltration: rainfall reduced by a saturation-excess
  //                    runoff term controlled by a_infil/b_infil, of which the
  //                    top layer accepts what accepted_inflow allows;
  //   * water_in_i>0 = what layer i accepts of the gravitational drainage of
  //                    the layer above, conductivity(theta_{i-1}, i-1), a
  //                    Clapp & Hornberger (1978) / Zeng & Decker (2009)
  //                    unsaturated hydraulic conductivity. A nearly full layer
  //                    accepts no more than it passes on, so the excess stays
  //                    above it and, at the surface, runs off;
  //   * water_out_i  = water_in_{i+1}, and for the deepest layer its own
  //                    gravitational drainage (water_flux[i] either way);
  //   * root_uptake_i= resource_depletion[i] (m yr^-1): the plants' uptake from
  //                    the layer, built from each strategy's
  //                    evapotranspiration_dt and divided by the patch area in
  //                    Patch::compute_rates.
  //
  // The final `aux_num` state slots accumulate diagnostic cumulative fluxes
  // (rainfall, infiltration, deep drainage, total root uptake, pulse runoff).
  //
  // This is an explicit, first-order representation; drainage is instantaneous
  // single-direction (no upward capillary flux between layers - that is handled
  // hydraulically inside the plant, by the leaf's root network, phylloptim's
  // E_from_Soil_to_Root_Collar).
  //
  // The inter-layer flux is FREE DRAINAGE, q = K(theta) evaluated at the donor
  // cell (Darcy's law with the matric gradient dropped), less what a nearly
  // full receiver refuses. Neither part contains an internode distance, so the
  // flux holds unchanged for unequal layers, and mass conservation,
  // sum(dz[i] * dtheta_i/dt) = infiltration - deep drainage - total uptake, is
  // exact for any widths (pinned on a graded profile in
  // test-tf24-water-budget.R) while no layer is held at the residual moisture.
  // A held layer's rate is set to zero, but layer_rates still passes its
  // drainage to the layer below and its uptake is still booked, so the budget
  // stops closing there.
  //
  // ⚠️ BUT IF ANYONE ADDS A CAPILLARY OR DIFFUSIVE FLUX -- and soil evaporation
  // will make them want to, because a top layer drying from above needs upward
  // flow -- THAT term does carry a distance, and it is NOT dz[i]. It is the
  // internode distance z_mid[i+1] - z_mid[i], and the conductivity at the
  // interface must be a distance-weighted harmonic or geometric mean of the two
  // cells, not either cell's K. Using dz[i] there would be correct only for equal
  // layers and would break the graded-profile budget test.
  void compute_rates(std::vector<S> const &resource_depletion)
  {
    using std::max;

    // Rainfall is floored at zero, as a counted backstop. Drivers are
    // interpolated with monotone slopes (Environment), so a read stays between
    // the values that bracket it and a non-negative series never reads
    // negative. Were it to, negative rainfall would give negative infiltration,
    // rain drying the soil, and at the residual moisture the guard below would
    // clamp the rate while `sum_rainfall` recorded the removal, so the water
    // budget would stop closing. The count says whether it ever binds.
    const double rainfall = rainfall_floored();
    const S infiltration =
      layer_rates(rainfall, vars.states, resource_depletion, vars.rates,
                  &water_flux, [this](int site) { note_clamp(site); });
    S total_resource_depletion = 0.0;
    for (size_t i = 0; i < soil_number_of_depths; i++) {
      resource_uptake[i] = resource_depletion[i];
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

  // The rainfall driver at `t`, by default the environment's time, floored at
  // zero for the reasons given on compute_rates(). A member, so that the soil's
  // sub-steps and TF24w's private root zones floor and count it as the
  // environment's own column does.
  double rainfall_floored() const { return rainfall_floored(time); }
  double rainfall_floored(double t) const {
    const double rainfall_raw = extrinsic_drivers.evaluate("rainfall", t);
    if (rainfall_raw < 0.0) {
      note_clamp(CLAMP_RAINFALL);
    }
    return std::max(0.0, rainfall_raw);
  }

  // The conductivity at theta held to [0, saturation], counting no clamp.
  template <class U> U conductivity(const U& theta, size_t layer) const {
    using std::max;
    using std::min;
    using std::pow;
    const double k_sat_layer = soil_parameter_value(K_sat_layers, K_sat, layer);
    const double soil_moist_sat_layer =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    const double n_psi_layer = soil_parameter_value(n_psi_layers, n_psi, layer);
    const U t = min(max(theta, U(0.0)), U(soil_moist_sat_layer));
    return k_sat_layer * pow(t / soil_moist_sat_layer, 2 * n_psi_layer + 3);
  }

  // The share of the rain the top layer takes at moisture theta, before it is
  // floored at zero: one less the saturation-excess runoff.
  template <class U> U infiltrated_share(const U& theta) const {
    using std::pow;
    const double soil_moist_sat_0 =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, 0);
    return U(1) - a_infil * pow(theta / soil_moist_sat_0, b_infil);
  }

  // The exponent b of the saturated share w in accepted_inflow, which sets how
  // near saturation a layer begins to refuse what it cannot pass on: w is
  // 1.5e-5 at half saturation and 1% at three-quarters. A flatter w refuses
  // water a layer has room for, and moves columns no layer of which is full;
  // a steeper one stiffens a layer held full above a slower one, where the
  // stiffness is about b (supply - out) / (theta_sat dz). 16 is where a
  // perched column's outcome stops depending on b. Hazard 12 has the
  // measurements.
  static constexpr double inflow_refusal_exponent = 16.0;
  // The most the saturated share w can reach, which it approaches only far
  // above saturation.
  static constexpr double inflow_refusal_bound = 100.0;
  // The sharpness p of smooth_min. The smooth minimum falls short of min(a, b)
  // most where a = b, by a factor 2^(-1/p), 16% at 4. A draining column's
  // interfaces sit near a = b, where its curvature grows with p, so a sharp
  // minimum puts a near-kink on the path of ordinary runs, which a difference
  // in a trait then straddles. The shortfall is weighted by the saturated
  // share, so it barely moves a column that passes on what it receives.
  static constexpr double inflow_limit_sharpness = 4.0;

  // A smooth minimum of a and b that never exceeds the smaller: the reciprocal
  // p-norm a b / (a^p + b^p)^(1/p), p = inflow_limit_sharpness. It is
  // symmetric, so it is evaluated from the larger argument, which keeps the
  // power in [0, 1] and finite however far apart the two are (conductivities
  // span ~25 orders of magnitude); the two orderings are one analytic
  // function, so choosing between them puts no kink in it. It is zero where
  // both are zero, its limit.
  template <class U> static U smooth_min(const U& a, const U& b) {
    using std::pow;
    const bool a_smaller = a < b;
    const U& lo = a_smaller ? a : b;
    const U& hi = a_smaller ? b : a;
    if (!(hi > 0.0)) {
      return U(0.0);
    }
    // The power is formed by squaring, so it stays real for a negative ratio
    // (an outflow an overfilled stage turned upward).
    static_assert(inflow_limit_sharpness == 4.0,
                  "smooth_min forms (lo / hi)^p by two squarings");
    const U r = lo / hi;
    const U r2 = r * r;
    return lo / pow(U(1.0) + r2 * r2, 1.0 / inflow_limit_sharpness);
  }

  // What layer `layer`, at moisture theta, takes of a supply from above when
  // it passes `out` on below. A dry layer takes the whole supply, into its free
  // pore space; a full one takes only what it passes on. Between the two, w is
  // read as the saturated share of the layer, which takes only what the layer
  // passes on while the rest takes freely:
  //
  //   accepted = supply - w (supply - smooth_min(supply, out)),
  //   w = W f / (W - 1 + f),  f = (theta / theta_sat)^b,
  //
  // with b = inflow_refusal_exponent and W = inflow_refusal_bound, so w is f
  // to within 1/W below saturation and exactly 1 at it.
  //
  // At saturation the layer's inflow is at most its outflow, so its rate is not
  // positive while its roots draw on it, and saturation is a barrier the
  // balance keeps by itself. Above saturation (an explicit stage can put it
  // there) w exceeds one and the layer refuses more, which pushes it back. What
  // a layer refuses stays in the layer above, which fills and refuses in turn;
  // at the surface it is runoff. Nothing is lost: each flux leaves one layer
  // and enters the next by the same amount.
  //
  // w is bounded above because, left as f, an explicit stage far above
  // saturation would weigh the refusal by theta^b, and the stages that follow
  // it would grow by that power until they overflowed, where a bounded flux
  // lets the step be rejected. The moisture is floored at zero in f only,
  // where f and its slope are already zero; a stage below zero is already
  // outside the domain.
  template <class U>
  U accepted_inflow(const U& supply, const U& out, const U& theta,
                    size_t layer) const {
    using std::max;
    using std::pow;
    const double sat =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    const U f = pow(max(theta, U(0.0)) / sat, inflow_refusal_exponent);
    constexpr double W = inflow_refusal_bound;
    const U w = W * f / (W - 1.0 + f);
    return supply - w * (supply - smooth_min(supply, out));
  }

  // Each layer's rate under `rainfall`, less the uptake u from it, held at the
  // residual moisture; returns the infiltration, the water that enters the top
  // layer. `note` counts each clamp, and `drainage`, if given, takes the flux
  // out of the base of each layer.
  //
  // A layer's inflow depends on what the layers below it pass on, so the
  // fluxes are found from the bottom up: the flux out of the deepest layer is
  // its free drainage, and the flux into each layer is accepted_inflow of what
  // the layer above offers (its free drainage, or the infiltrated rain at the
  // surface) against what this layer passes on.
  template <class U, class Note>
  U layer_rates(double rainfall, const std::vector<U>& theta,
                const std::vector<U>& u, std::vector<U>& rate,
                std::vector<U>* drainage, Note note) const {
    using std::max;
    const size_t n = soil_number_of_depths;
    const U share = infiltrated_share(theta[0]);
    // Where the top layer is wet enough that this goes negative, the infiltrated
    // fraction stops reading that layer's state at all.
    if (share < U(0.0)) {
      note(CLAMP_INFILTRATION);
    }
    // An explicit stage can evaluate theta outside [0, saturation], where the
    // conductivity curve is held.
    auto free_drainage = [&](size_t i) {
      const U& t = theta[i];
      if (t < U(0.0) ||
          t > U(soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, i))) {
        note(CLAMP_SOIL_CONDUCTIVITY);
      }
      return conductivity(t, i);
    };
    U out = free_drainage(n - 1);
    for (size_t i = n; i-- > 0;) {
      const U supply =
        i == 0 ? U(rainfall * max(U(0.0), share)) : free_drainage(i - 1);
      const U in = accepted_inflow(supply, out, theta[i], i);
      if (drainage != nullptr) {
        (*drainage)[i] = out;
      }
      rate[i] = (in - out - u[i]) / dz[i];
      // At or below the residual moisture a layer is not dried further, and a
      // NaN rate (from the retention curve near residual) is held too.
      if (theta[i] <= soil_moist_residual && !(rate[i] > 0.0)) {
        note(CLAMP_SOIL_POSITIVITY);
        rate[i] = U(0.0);
      }
      out = in;
    }
    return out;
  }

  // The times in (from, to) where the rain changes form: the soil's knots.
  std::vector<double> rainfall_knots(double from, double to) const {
    return extrinsic_drivers.knots("rainfall", from, to);
  }

  // The layers' rates at `time_` under a given uptake u from each, as
  // compute_rates finds them.
  template <class U>
  void layer_rates(double time_, const std::vector<U>& theta,
                   const std::vector<U>& u, std::vector<U>& rate) const {
    rate.resize(theta.size());
    layer_rates<U>(rainfall_floored(time_), theta, u, rate, nullptr,
                   [this](int site) { note_clamp(site); });
  }


  // convert soil moisture to soil water potential
  S psi_from_soil_moist(S soil_moist_, size_t layer) const {
    return psi_from_soil_moist_noting<true>(soil_moist_, layer);
  }

  // The same potential with no clamp counted, for a record of the state: reading
  // the state must not change what the run reports it clamped.
  S psi_from_soil_moist_unnoted(S soil_moist_, size_t layer) const {
    return psi_from_soil_moist_noting<false>(soil_moist_, layer);
  }

  template <bool Note>
  S psi_from_soil_moist_noting(S soil_moist_, size_t layer) const {
    using std::max;
    using std::min;
    using std::pow;
    // Floor at the residual moisture: the retention curve (negative exponent)
    // diverges to +inf as theta->0, so an empty layer would otherwise yield a
    // non-finite potential. At/below theta_r the potential is large but finite
    // and the plant's root vulnerability curve has already shut uptake to ~0.
    if (Note && soil_moist_ < S(soil_moist_residual)) {
      note_clamp(CLAMP_SOIL_MOISTURE_FLOOR);
    }
    const S t = max(soil_moist_, S(soil_moist_residual));
    const double a_psi_layer = soil_parameter_value(a_psi_layers, a_psi, layer);
    const double soil_moist_sat_layer =
      soil_parameter_value(soil_moist_sat_layers, soil_moist_sat, layer);
    const double n_psi_layer = soil_parameter_value(n_psi_layers, n_psi, layer);
    const S psi =
      a_psi_layer * pow(t / soil_moist_sat_layer, -n_psi_layer);
    // Cap at a large-but-finite ceiling: near theta_r the retention curve
    // gives psi ~ 1e8 MPa, which pushes the leaf hydraulic solve non-finite and
    // corrupts the soil feedback. Root conductivity is already negligible at the
    // ceiling (see soil_psi_max_), so the cap keeps every downstream
    // calculation finite without changing uptake.
    if (Note && psi > S(soil_psi_max_)) {
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
    return pow((psi / a_psi_layer), (-1 / n_psi_layer)) * soil_moist_sat_layer;
  }

  double soil_moist_from_psi(double psi) const {
    return soil_moist_from_psi(psi, 0);
  }

  // Each is read once per individual per ODE derivs evaluation, always at the
  // current `time`, so the memo collapses the repeated string-keyed lookups
  // across individuals into one lookup per driver per distinct time. Every one
  // returns what a direct `extrinsic_drivers_evaluate(name, time)` returns.
  // Incoming shortwave at the surface, horizontal plane, above the canopy (W m^-2)
  double get_SWdown()     const { return driver_at(DRIVER_SWDOWN); }
  double get_atm_vpd()    const { return driver_at(DRIVER_ATM_VPD); }
  double get_ca()         const { return driver_at(DRIVER_CA); }
  double get_leaf_temp()  const { return driver_at(DRIVER_LEAF_TEMP); }
  double get_atm_o2_kpa() const { return driver_at(DRIVER_ATM_O2_KPA); }
  double get_atm_kpa()    const { return driver_at(DRIVER_ATM_KPA); }
  // Above-canopy wind speed U0 (m s^-1), Penman-Monteith aerodynamic resistance
  double get_wind_speed() const { return driver_at(DRIVER_WIND_SPEED); }


  std::vector<double> get_soil_water_state() const {
    std::vector<double> out = rebind(vars.states);
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
    const std::vector<double> all = rebind(vars.states);
    return {all.end() - aux_num, all.end()};
  }
  // Cumulative depth to the BOTTOM of each layer, m.
  std::vector<double> get_soil_depths() const { return z; }
  // The WIDTH of each layer, m -- what set_soil_layer_widths() takes, and what a
  // caller needs to turn the per-layer theta into a water volume: storage is
  // sum(theta * dz), never sum(theta) * dz[1].
  std::vector<double> get_soil_layer_widths() const { return dz; }


  Internals<double> r_internals() const {
    Internals<double> out(vars.state_size(), vars.aux_size(), vars.resource_size());
    out.states = rebind(vars.states);
    out.rates = rebind(vars.rates);
    return out;
  }

  // The generic resource-pulse hook: TF24's resources are its soil
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
    //
    // ⚠️ SO IS WHAT THE LAYER TAKES. The branch is chosen on passive values,
    // but a pulse that fills the layer takes the room it had, which is a
    // function of the state: the layer ends at saturation whatever it held, so
    // its row in the state it replaced is zero, and infiltration and runoff
    // move with that state, by -dz and +dz. A pulse taken whole, or shed whole
    // by a layer already full, adds a constant. Taken passive, the filled
    // branch would carry the identity's rows, and the stand gradient across
    // the pulse would be finite and wrong.
    const S theta = vars.state(layer);
    const S room = pulse_room(theta, layer);
    const double capacity = std::max(0.0, odelia::util::to_passive(room));
    const S accepted = !(capacity < depth) ? S(depth)
                       : capacity > 0.0    ? room
                                           : S(0.0);
    const S excess = S(depth) - accepted;

    vars.set_state(layer, theta + accepted / dz[layer]);
    // Direct state increments, not rates: the trailing slots are integrated
    // from their rates during a leg, and a pulse happens between legs.
    const size_t n = soil_number_of_depths;
    vars.set_state(n,     vars.state(n)     + S(depth));     // sum_rainfall
    vars.set_state(n + 1, vars.state(n + 1) + accepted);     // sum_infiltration
    vars.set_state(n + 4, vars.state(n + 4) + excess);       // sum_pulse_runoff

    psi_soil_valid_ = false;
    return {odelia::util::to_passive(accepted),
            odelia::util::to_passive(excess)};
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
    initial_states = rebind(vars.states);
    psi_soil_valid_ = false;
}

  // What a run reads from this environment, as plain values a saved scenario
  // can keep: the soil column, its parameters, the water a run starts from and
  // the drivers. A run resets the time, the light field and the soil water to
  // the start, so those are not kept, and the water kept is the starting water.
  // set_recipe() is the inverse: recipe() of a fresh environment it sets from a
  // complete recipe is identical() to that recipe.
  Rcpp::List recipe() const {
    SEXP layered = R_NilValue;
    if (use_layered_soil_parameters) {
      layered = Rcpp::List::create(_["soil_moist_sat"] = soil_moist_sat_layers,
                                   _["K_sat"] = K_sat_layers,
                                   _["a_psi"] = a_psi_layers,
                                   _["n_psi"] = n_psi_layers);
    }
    const std::vector<double> water(initial_states.begin(),
                                    initial_states.begin() + soil_number_of_depths);
    return Rcpp::List::create(_["layer_widths"] = dz,
                              _["layer_depths"] = z,
                              _["depth"] = depth,
                              _["soil_moist_sat"] = soil_moist_sat,
                              _["K_sat"] = K_sat,
                              _["a_psi"] = a_psi,
                              _["n_psi"] = n_psi,
                              _["layered_soil_parameters"] = layered,
                              _["a_infil"] = a_infil,
                              _["b_infil"] = b_infil,
                              _["f_PAR"] = f_PAR,
                              _["water_mixing_rate"] = water_mixing_rate,
                              _["soil_water"] = water,
                              _["drivers"] = r_driver_recipe());
  }

  // The column is set from its widths, which checks them, and then given the
  // recipe's own layer depths and total, which may differ from the running sum
  // of the widths in the last bits: set_soil_number_of_depths() builds them
  // from the total, and TF24's numbers move with that difference.
  //
  // Made for a fresh environment, which is the only one plant hands it. There a
  // field the recipe lacks takes its default, so a recipe saved before a field
  // was added still sets an environment. A recipe that sets the column resets
  // what belongs to the old one, as the geometry setters do: the water starts
  // at half the saturation in force, as the constructor's does, and any
  // per-layer parameters the recipe does not give again are dropped. A name the
  // recipe has and this does not know is refused, so a misspelt field is not
  // dropped unseen.
  void set_recipe(Rcpp::List r) {
    static const std::vector<std::string> known = {
      "layer_widths", "layer_depths", "depth", "soil_moist_sat", "K_sat",
      "a_psi", "n_psi", "layered_soil_parameters", "a_infil", "b_infil",
      "f_PAR", "water_mixing_rate", "soil_water", "drivers"};
    if (r.size() > 0 && Rf_isNull(r.names())) {
      util::stop("A TF24 environment's recipe must be a named list");
    }
    const Rcpp::CharacterVector given =
      r.size() == 0 ? Rcpp::CharacterVector() : Rcpp::CharacterVector(r.names());
    for (R_xlen_t i = 0; i < given.size(); ++i) {
      const std::string name = Rcpp::as<std::string>(given[i]);
      if (std::find(known.begin(), known.end(), name) == known.end()) {
        util::stop("A TF24 environment's recipe has no field '" + name + "'");
      }
    }
    auto has = [&r](const char* name) { return r.containsElementNamed(name); };
    auto set_double = [&](const char* name, double& field) {
      if (has(name)) {
        field = Rcpp::as<double>(r[name]);
      }
    };

    const bool column = has("layer_widths");
    if (!column && (has("layer_depths") || has("depth"))) {
      util::stop("A TF24 environment's recipe gives layer_depths or depth "
                 "without the layer_widths they belong to");
    }
    if (column) {
      set_soil_layer_widths(Rcpp::as<std::vector<double>>(r["layer_widths"]));
      const std::vector<double> depths = has("layer_depths")
        ? Rcpp::as<std::vector<double>>(r["layer_depths"]) : z;
      const double total = has("depth") ? Rcpp::as<double>(r["depth"]) : z.back();
      if (depths.size() != z.size()) {
        util::stop("The recipe's layer_depths has " +
                   util::to_string(depths.size()) + " entries for " +
                   util::to_string(z.size()) + " layer_widths");
      }
      // A running sum of n widths drifts from the product built from the total
      // by about one rounding per layer, so the tolerance grows with n.
      const double tol = 4.0 * static_cast<double>(z.size() + 1) *
        std::numeric_limits<double>::epsilon();
      auto close = [tol](double a, double b) {
        return std::abs(a - b) <= tol * std::max(std::abs(a), std::abs(b));
      };
      for (size_t i = 0; i < z.size(); ++i) {
        if (!close(depths[i], z[i])) {
          util::stop("The recipe's layer_depths are not the running sum of its "
                     "layer_widths: layer " + util::to_string(i + 1) +
                     " ends " + util::format_double(depths[i] - z[i]) +
                     " m from the sum of the widths above it");
        }
      }
      if (!close(total, z.back())) {
        util::stop("The recipe's depth is " +
                   util::format_double(total - z.back()) +
                   " m from the bottom of its last layer");
      }
      z = depths;
      depth = total;
      finalise_soil_geometry_();
    }

    set_double("soil_moist_sat", soil_moist_sat);
    set_double("K_sat", K_sat);
    set_double("a_psi", a_psi);
    set_double("n_psi", n_psi);
    set_double("a_infil", a_infil);
    set_double("b_infil", b_infil);
    set_double("f_PAR", f_PAR);
    set_double("water_mixing_rate", water_mixing_rate);
    if (has("layered_soil_parameters")) {
      SEXP layered = r["layered_soil_parameters"];
      if (Rf_isNull(layered)) {
        use_layered_soil_parameters = false;
        soil_moist_sat_layers.clear();
        K_sat_layers.clear();
        a_psi_layers.clear();
        n_psi_layers.clear();
      } else {
        // A layer value left out takes the scalar field for every layer, as
        // set_tf24_soil() reads NULL.
        const Rcpp::List l(layered);
        auto layer = [&l](const char* name) -> SEXP {
          return l.containsElementNamed(name) ? SEXP(l[name]) : R_NilValue;
        };
        set_soil_parameters(soil_number_of_depths, layer("soil_moist_sat"),
                            layer("K_sat"), layer("a_psi"), layer("n_psi"));
      }
    }
    if (has("soil_water")) {
      set_soil_water_state(Rcpp::as<std::vector<double>>(r["soil_water"]));
    } else if (column) {
      set_soil_water_state(
        std::vector<double>(soil_number_of_depths, soil_moist_sat * 0.5));
    }
    if (has("drivers")) {
      r_set_driver_recipe(r["drivers"]);
    }
    // The memo is keyed by the time alone, which this leaves where it was.
    drivers_ = driver_memo{};
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
  //
  // A layer set above its saturation is refused here, where every run takes
  // its starting soil, rather than when set: the saturation can be changed
  // after the moisture is, and a run could not leave such a state anyway.
  virtual void clear_state() {
    util::check_length(initial_states.size(), vars.state_size());
    const std::string refusal = saturation_refusal(initial_states.begin());
    if (!refusal.empty()) {
      util::stop("The run cannot start: " + refusal +
                 ". Set the soil water state at or below saturation");
    }
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

    // Each layer's water potential from its own retention curve, floor and cap.
    Rcpp::NumericVector rcpp_soil_psi_vec(soil_moist_list.size());
    for (size_t i = 0; i < soil_moist_list.size(); ++i) {
      rcpp_soil_psi_vec[i] = odelia::util::to_passive(
        psi_from_soil_moist_unnoted(S(soil_moist_list[i]), i));
    }

    auto const &soil_moist_cumulative_flux_list = get_soil_water_state_cumulative_flux();
    auto rcpp_soil_moist_vec_cumulative_flux = Rcpp::NumericVector(soil_moist_cumulative_flux_list.begin(), soil_moist_cumulative_flux_list.end());

    // Layer WIDTHS as well as boundaries. Reported because `soil_moist` is theta,
    // an intensive quantity: turning it into a water volume needs each layer's
    // width, not soil_depth[[1]], which is right only for equal layers.
    auto const &soil_width_list = get_soil_layer_widths();
    auto rcpp_soil_width_vec = Rcpp::NumericVector(soil_width_list.begin(), soil_width_list.end());

    return Rcpp::List::create(
        // auto ret = get_state(environment.extrinsic_drivers, time);

        _["light_availability"] = light_availability.r_get_state(),
        _["soil_moist"] = rcpp_soil_moist_vec,
        _["soil_psi"] = rcpp_soil_psi_vec,
        _["soil_depth"] = rcpp_soil_depth_vec,
        _["soil_layer_width"] = rcpp_soil_width_vec,
        _["soil_moist_cumulative_flux"] = rcpp_soil_moist_vec_cumulative_flux
    );
  }
  };
}

#endif
