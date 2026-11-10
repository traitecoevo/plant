// -*-c++-*-
// TF24_Strategy's leaf coupling: net production through the leaf's solve,
// the leaf's water uptake, and the rows the gradient path records for it.
// The class is declared in tf24_model.h, which includes this file.
#ifndef PLANT_PLANT_TF24_LEAF_H_
#define PLANT_PLANT_TF24_LEAF_H_

#include <plant/models/tf24_model.h>

namespace plant {

// Reads the soil consumption the leaf's last solve left, so the leaf must have
// been solved for this plant at this step.
template <typename S>
S TF24_Strategy<S>::evapotranspiration_dt(const S& area_leaf_, int soil_layer) {
  if constexpr (std::is_same_v<S, double>) {
    return leaf.soil_consumption_[soil_layer] * area_leaf_;
  } else {
    return leaf_soil_consumption_[soil_layer] * area_leaf_;
  }
}

// Two of the leaf's outputs re-enter the active chain, and they reach the
// operating point differently.
//
// PROFIT is the objective at the point the leaf chose, so by the envelope
// theorem its response at an interior optimum is the direct one and the point's
// own movement drops out. PER-LAYER UPTAKE is set as a side effect AT that
// point, so it consumes the choice rather than being it and the movement is part
// of its answer.
//
// The leaf answers at whatever scalar it is asked at, so this hands it the
// inputs and takes the outputs back on the tape. Only one derivative crosses as
// a number -- the gradient of the condition that places an interior collar --
// and it crosses because reverse mode cannot nest a tangent above its own
// scalar, where forward mode can.
template <typename S>
void TF24_Strategy<S>::record_leaf_outputs(const S& PPFD_int,
                                           const std::vector<S>& psi_soil,
                                           const S& conductance_max,
                                           double NIR_per_PPFD_int) {
  using phylloptim::Leaf;
  const std::size_t n_layer = psi_soil.size();

  // The two _25 traits and dark respiration enter as the temperature-adjusted
  // values the leaf's kernels take. The Arrhenius is linear in its reference
  // value, so the chain from the trait is the factor between them and the tape
  // carries the rest.
  // ⚠️ THE TRAITS THEMSELVES, AND EVERY SLOT NAMED. The pack is the leaf model's
  // own parameter enumeration, so the temperature responses, both Weibull scale
  // parameters and both critical potentials are DERIVED inside the leaf from
  // these -- which is what puts their chains on this tape instead of asking this
  // caller to hand over four numbers it would have to keep consistent.
  //
  // Filled by NAMED INDEX rather than as an aggregate: an initialiser shorter
  // than the array zero-fills the rest with no diagnostic at all, and this list
  // has been the wrong length twice.
  //
  // ⚠️ SEATED FIRST, THEN OVERRIDDEN. Starting from the leaf's OWN pack means no
  // slot is left at a value this caller did not choose: the parameters of the
  // seven cost curves this strategy does not run are carried as the leaf's
  // constants rather than as zeros, which is what a bare fill() would have made
  // them -- a wrong VALUE, not merely a missing row, in any expression that
  // reads one.
  phylloptim::leaf_pars<S> in;
  const phylloptim::leaf_pars<double> seated = leaf.passive_pars();
  for (std::size_t i = 0; i < in.size(); ++i) {
    in[i] = S(seated[i]);
  }
  // What THIS strategy owns, and therefore what carries rows.
  in[phylloptim::par_vcmax_25] = pars.vcmax_25;
  in[phylloptim::par_stem_c] = pars.stem_c;
  in[phylloptim::par_stem_P50] = pars.stem_P50;
  in[phylloptim::par_root_c] = pars.root_c;
  in[phylloptim::par_root_P50] = pars.root_P50;
  in[phylloptim::par_TF24_beta2] = pars.TF24_beta2;
  in[phylloptim::par_jmax_25] = pars.jmax_25;
  in[phylloptim::par_phi_J] = pars.phi_J;
  in[phylloptim::par_curv_fact_elec_trans] = pars.curv_fact_elec_trans;
  in[phylloptim::par_curv_fact_colim] = pars.curv_fact_colim;
  in[phylloptim::par_TF24_cost_scale] = pars.TF24_cost_scale;
  in[phylloptim::par_R_d_25] = pars.R_d_25;
  in[phylloptim::par_kmax] = conductance_max;
  in[phylloptim::par_TF24_floor_lambda_o] = pars.TF24_floor_lambda_o;
  // The intercepted PPFD this cohort stands in, and the one route every other cohort's
  // height and leaf area have into its carbon: the electron transport reads
  // this slot. Its VALUE is already seated -- optimise_at strips the same
  // number to passive for the leaf -- so what this line adds is the row, on
  // exactly the footing `kmax` is on.
  in[phylloptim::par_PPFD] = PPFD_int;
  // The near-infrared is a fixed multiple of the same light (optimise_at seats
  // its value), and with the energy balance on it heats the leaf: its row is
  // light's other route to leaf temperature.
  in[phylloptim::par_NIR] = NIR_per_PPFD_int * PPFD_int;

  std::vector<S> r_R_H_min, r_R_V_sum;
  // Carbon becomes resistance HERE, at this scalar, because the architecture
  // model is this strategy's -- the leaf takes the resistances and knows nothing
  // about carbon. So the chain lands on the same tape as everything else rather
  // than crossing as a hand-written d(uptake)/d(carbon).
  // The thickness is read off the leaf's own soil profile rather than derived
  // again here: the vertical resistance scales with dz^2, so two definitions
  // drifting apart would be a silent squared factor neither side could see.
  // ⚠️ THROUGH layer_thicknesses, THE ONE DEFINITION OF dz, and off the LEAF's
  // own profile rather than a second copy of the environment's. The vertical
  // resistance scales with dz^2, so two definitions drifting apart would be a
  // silent SQUARED factor that neither side could detect.
  phylloptim::root_network_from_carbon<S>(
      root_carbon_per_leaf_area_,
      phylloptim::layer_thicknesses(leaf.roots_.soil_depth_),
      beta_R_H, beta_R_V, r_R_H_min, r_R_V_sum);
  // The supply as a VIEW over what this caller owns. The curve arrives as
  // (root_P50, root_c) -- root_b is derived and a row in it reaches nothing.
  const phylloptim::SupplyInputs<S> supply{psi_soil, r_R_H_min, r_R_V_sum,
                                       in[phylloptim::par_root_P50],
                                       in[phylloptim::par_root_c]};

  // The leaf clamps in double on both paths, so which path a clamp fired on is a
  // question of when rather than of the scalar, and a delta answers it. Taken on
  // the throwing exit too, because a refusal is still a pass over the leaf.
  const leaf_clamp_tally clamps_before = leaf_clamps();
  // ⚠️ THE SAME NUMBER THE INTERIOR CLOSURE DIVIDES BY, so a guard that tests a
  // different one passes while the closure divides by a bad slope.
  // profit_curvature_at_collar() is dM/dp as a first derivative of the marginal the
  // solve roots, and the rows are taped from M inside operating_collar.
  // The marginal at the collar the solve left, for a refusal message. On a COPY
  // because dprofit_at_collar_psi drives the model to that collar and this sits
  // on a path that is about to throw with the leaf's own state still readable.
  auto marginal_here = [&]() -> double {
    phylloptim::Leaf probe = leaf;
    bool ok = false;
    const double m = probe.dprofit_at_collar_psi<
        phylloptim::Leaf::CostCurve::TF24_floor>(leaf.opt_root_psi_, &ok);
    return ok ? m : std::numeric_limits<double>::quiet_NaN();
  };
  double slope = std::numeric_limits<double>::quiet_NaN();
  Leaf::BilledOutputs<S> got;
  try {
    if (leaf.operating_point_kind() == Leaf::OperatingPointKind::Interior) {
      slope = leaf.template profit_curvature_at_collar<
          phylloptim::Leaf::CostCurve::TF24_floor>();
      // The interior derivation divides by the profit's curvature. At a bound
      // the slope is the bound's own and this floor is not about it.
      note_curvature(slope);

      // profit_curvature_at_collar is a difference of the marginal, measured
      // rather than assembled analytically (phylloptim says why), so there is no
      // independent curvature to check it against here.
      //
      // Two findings, two messages: DO NOT merge them into one sentence. They
      // are not the same KIND of problem -- the floor limb is a real degeneracy
      // of the model, the sign limb is a SOLVER failure -- and one sentence
      // reports a positive curvature of 34.4 "against a floor of 0.001", four
      // orders of magnitude ABOVE that floor, which reads as a magnitude
      // failure and is not one.
      if (!(slope < 0.0)) {
        // ⚠️ A POSITIVE CURVATURE HERE IS EVIDENCE AGAINST THE CURVATURE, not
        // against the operating point. An interior collar is reached only on a
        // bracket the marginal crosses downward, so the curvature at it must be
        // non-positive -- and the collar solve now probes dry of every root it
        // returns and counts the exceptions. On the 105-year stand in
        // scripts/profile-stand-gradient.R it found the marginal rising ZERO times
        // in 2,829,445 interior solves, while this reported +34.4; a centred
        // difference of the solve's own marginal at that collar gave -9.63.
        //
        // So the two spellings of the marginal disagree, and the one that matches
        // the difference is the solve's. `condition.marginal` and
        // `leaf.collar_resid_` are those two spellings, carried here because their
        // disagreement is the finding -- both must be ~0 at an interior point, and
        // where they are not, profit_at is not the function that placed this
        // collar and its second derivative belongs to something else.
        //
        // format_double, not to_string: these numbers run to ~1e-6 and below,
        // where to_string's six DECIMAL places print "0.000000". format_double
        // keeps six SIGNIFICANT figures and exists in util.h for this.
        refuse(std::string("TF24 gradient: the profit curvature at this "
                           "operating point is ") +
               util::format_double(slope) +
               ", which is not negative -- so the interior derivation cannot "
               "divide by it. The marginal at this collar is " +
               util::format_double(marginal_here()) +
               "; the collar is " +
               util::format_double(leaf.opt_root_psi_) +
               ". An interior collar sits on a downward crossing, so a positive "
               "curvature here indicts this derivation rather than the point. See "
               "docs/design/leaf-derivatives.md");
      } else if (std::abs(slope) < curvature_floor()) {
        // format_double for the same reason: a curvature inside a floor of 1e-3 is
        // where to_string's six decimal places start losing the number.
        refuse("TF24 gradient: the leaf's profit curvature at this operating "
               "point is " + util::format_double(slope) +
               ", inside the floor of " + util::format_double(curvature_floor()) +
               ", so the interior derivation would divide by too little to "
               "answer");
      }
    }
    // Placed, then evaluated. Which condition placed it is the leaf's to say;
    // whether the curvature it divides by is usable was this caller's, above.
    // The supply, taken ONCE at the point the solve left. The soil reaches
    // everything the leaf answers with through this and nothing else, so both
    // calls below read it instead of re-recording the same Ohm's law over the
    // same tabulated integral.
    const phylloptim::Leaf::SupplyDraw<S> draw =
        leaf.template supply_draw_at<S>(S(leaf.opt_root_psi_), supply);
    // `slope` is NaN unless the point is interior, which is the one kind that
    // divides by it -- so this hands over the number the guard above refused on
    // and leaves every other kind to close on its own condition.
    // `held` carries the stem potential and ci the interior arm solved on its way
    // to the collar, so billing reads them rather than solving both again.
    Leaf::HeldCoords<S> held;
    const S collar =
        leaf.template operating_collar<phylloptim::Leaf::CostCurve::TF24_floor, S>(
            draw, in, slope, &held);
    got = leaf.template billed_outputs_at<phylloptim::Leaf::CostCurve::TF24_floor,
                                          S>(collar, draw, in, &held);
  } catch (const std::runtime_error& e) {
    note_leaf_clamps(clamps_before);
    // The leaf produced nothing to record a row against, so every output takes
    // the value its own solve fixed and carries none -- which is what a refusal
    // costs everywhere else here, and lets the sweep reach the poll that reads
    // this.
    refuse(std::string("TF24 gradient: ") + e.what());
    leaf_carbon_ = S(leaf.carbon());
    leaf_soil_consumption_.assign(n_layer, S(0.0));
    for (std::size_t j = 0; j < n_layer; ++j) {
      leaf_soil_consumption_[j] = S(leaf.soil_consumption_[j]);
    }
    return;
  }
  note_leaf_clamps(clamps_before);

  // A collar the theorem could not place costs every output that reads it, and
  // profit is not one of them: at an interior point it reads the collar fixed. An
  // unplaceable collar throws, and the catch above refuses at that grain.

  // ⚠️ THE VALUE IS THE LEAF'S, NOT THIS COMPOSITION'S. The solve is what plant's
  // rates read, so a second copy of it could only disagree; what is taken here is
  // the derivative, recorded against a value the forward pass already fixed.
  //
  // The objective is profit, and it is the one output a refusal does not stop
  // being asked for: it is built by the envelope theorem and reads no curvature,
  // so it still carries a row at a point where the water rows have none.
  auto carry = [&](double value, const S& from, S& into, bool objective) -> void {
    if (recorded_refusal->happened() && !objective) {
      into = S(value);
      return;
    }
    // Named rather than braced: the rows cross as a span, which borrows storage
    // instead of owning it, so a stack array is what holds them here -- and is
    // what the span was for, since a braced list built a vector per call.
    const odelia::input_and_derivative<S> rows[]{{from, 1.0}};
    const odelia::supply_report report =
        odelia::supply_derivatives<S>(value, rows, into);
    if (report.whole) {
      return;
    }
    refuse(std::string("TF24 gradient: `") + (objective ? "profit" : "uptake") +
           "` cannot be recorded at a point the solve called " +
           Leaf::operating_point_kind_name(leaf.operating_point_kind()) + ": " +
           report.why);
  };

  // The carbon kept is carried where the objective was, first and as the
  // objective is: at the default price its rows ARE the objective's envelope
  // rows (the price adds an exact zero), and the refusals that would cost it the
  // collar's movement throw before here, into the value-only branch above.
  carry(leaf.carbon(), got.carbon, leaf_carbon_, true);
  leaf_soil_consumption_.assign(n_layer, S(0.0));
  util::check_length(got.uptake.size(), n_layer);
  for (std::size_t j = 0; j < n_layer; ++j) {
    carry(leaf.soil_consumption_[j], got.uptake[j], leaf_soil_consumption_[j],
          false);
  }
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
template <typename S>
S TF24_Strategy<S>::net_mass_production_dt(const TF24_Environment<S>& environment,
                                const S& height, const S& area_leaf_,
                                const S& height_inverse,
                                const S& sapwood_departure,
                                const S& leaf_departure,
                                const S& tracked_marginal) {
  // height_inverse (= 1/height) is supplied by the shared individual.h interface
  // (cached aux); unused here as the TF24 root-water path works in height directly.
  (void)height_inverse;
  const S mass_leaf_    = mass_leaf(area_leaf_);
  const S area_sapwood_ = area_sapwood(area_leaf_, sapwood_departure);
  const S mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const S area_bark_    = area_bark(area_leaf_);
  const S mass_bark_    = mass_bark(area_bark_, height);
  const S mass_root_    = mass_root(area_leaf_);

  int soil_number_of_depths_ = environment.get_soil_number_of_depths();
  const std::vector<double>& soil_depths_ = environment.z;



  // The intercepted light that drives the leaf optimisation depends on the shading
  // model and is computed below (just before the optimisation), once the
  // depth-independent inputs are ready.

  // psi_soil (-MPa): the patch column's, computed once per soil state and cached
  // in the environment -- or, in TF24w, this cohort's own column (#681).
  const std::vector<S>& psi_soil = soil_water_potential_for_leaf(environment);
  
// find leaf specific max hydraulic conductance (kg m^-2 LA s^-1 MPa ^-1)
  // pars.K_s: sapwood-specific conductivity of the TERMINAL segment -- which is
  //   the whole-stem value while pars.D_c == 0 (kg m^-1 s^-1 MPa^-1)
  // pars.theta: huber value
  // eta_c: accounts for average position of leaf mass
  // height: maximum plant height
  //
  // The flow path runs from the base to the leaf-area-weighted mean leaf
  // height, height*eta_c, NOT to the apex; eta_c therefore scales the UPPER
  // LIMIT of the path integral rather than the resistance. The two readings
  // differ by the constant eta_c^-beta, which is H-independent and so degenerate
  // with K_s -- nothing downstream can distinguish them, and the K_s
  // reparameterisation absorbs the difference entirely.
  //
  // Setting D_c, theta_c and L_tip all to zero makes effective_path_length
  // return height*eta_c having performed no arithmetic, so this expression is
  // then bit-identical to the height-linear code it replaces.
  const S stem_path_length = stem_hydraulics::effective_path_length(
      S(height * eta_c), pars.L_tip, stem_path_exponent());
  // ⚠️ NUMERATOR IS THE STATE, NOT pars.theta (#516). Supply per unit leaf area
  // follows the Huber value the plant actually carries, so a stem that keeps its
  // conducting area while the canopy thins supplies each remaining leaf better.
  // That is the acclimation half of #516, and it is also what the sapwood
  // controller differentiates: psi enters here as exp(psi), which is why
  // d(kmax)/d(psi) is just kmax. Stored so the controller reads the SAME value
  // rather than recomputing it -- recomputing it is how the two silently drift
  // apart when the denominator changes, as it did in #617.
  const S leaf_specific_conductance_max =
    pars.K_s * sapwood_per_leaf_area(sapwood_departure) / stem_path_length;
  leaf_specific_conductance_max_ = leaf_specific_conductance_max;

  // Fine-root carbon is distributed over depth using the same cumulative shape
  // function Q() used for the leaf canopy, but parameterised over soil depth
  // instead of crown height. Q(z, rooting_depth, 0.2) gives the fraction of
  // roots *below* depth z, so layer a holds root_mass_scale * (Q(z_{a-1}) -
  // Q(z_a)). rooting_depth is capped at the soil column depth and the loop
  // breaks once Q reaches 0, so deeper layers stay at zero.
  //
  // The quantity is carbon PER UNIT LEAF AREA, and that is a hard requirement
  // rather than a convenience: the leaf is intensive, it takes resistances built
  // from this, and a network built from absolute carbon is five vectors of
  // positive numbers that no check on the far side can distinguish -- uptake
  // then comes back wrong by the leaf area with nothing raised. mass_root() is
  // strictly linear in area_leaf, so dividing it back out is exactly
  // root_mass_carbon_scale * pars.a_r1 and costs no arithmetic.
  //
  // Reuse the member buffer (assign refills + zeroes without reallocating when
  // the layer count is unchanged); zeroing matters because the loop below breaks
  // early below the rooting depth, leaving deep layers that must read as 0.
  // TODO (perf): the scale (83.26) is hard-coded and should become a trait.
  root_carbon_per_leaf_area_.assign(soil_number_of_depths_, 0.0);
  root_carbon_value_.assign(soil_number_of_depths_, 0.0);

  // Counted, and this one severs a state row rather than smoothing a number:
  // above the cap the whole root profile stops depending on height, so every
  // layer's carbon -- and the uptake built on it -- reads a constant.
  if (height > pars.rooting_depth_max) {
    note_clamp(CLAMP_ROOTING_DEPTH);
  }
  const S rooting_depth = std::min(height, pars.rooting_depth_max);
  const S root_mass_scale = root_mass_carbon_scale * pars.a_r1;

  {
    using odelia::util::to_passive;
    S prev_q = 1.0;
    for (int a = 0; a < soil_number_of_depths_; ++a) {
      if (prev_q == 0) {
        break;
      }
      const S q = Q(soil_depths_[a], rooting_depth, pars.root_depth_shape_eta);
      root_carbon_per_leaf_area_[a] = root_mass_scale * (prev_q - q);
      root_carbon_value_[a] = to_passive(root_carbon_per_leaf_area_[a]);
      prev_q = q;
    }
  }

  // Carbon -> resistance: the root architecture model runs here, on the plant
  // side. The leaf's supply solve reads r_R_H_min and r_R_V_sum and knows
  // nothing about root carbon, the horizontal/vertical split or the layer
  // thickness, so this strategy owns the model exactly as it already owns the
  // conductance-versus-height model above, and hands over the reduced quantity.
  // The widths are phylloptim's layer_thicknesses of z, which the environment
  // computes once when its geometry is set (z_thickness) -- do not open-code
  // them, since the vertical resistance scales with dz^2 and the two sides
  // drifting apart would be a silent squared factor neither side could detect.
  //
  // TODO(next rebaseline): pass environment.dz, the widths the environment was
  // given, instead of differencing `z`. The two agree only to rounding (on the
  // default 1.5 m over 5, layers 3-5 differ in the last bit), so the switch
  // moves every TF24 number.
  phylloptim::root_network_from_carbon(
      root_carbon_value_, environment.z_thickness, beta_R_H,
      beta_R_V, root_network_);

  // Reuse geometry precomputed by environment; avoids rebuilding z midpoints each call.
  leaf.roots_.z_soil_mid_ = environment.get_soil_mid_depths();
  leaf.roots_.use_precomputed_z_soil_mid_ = true;

  // Per-timestep above-canopy wind, read only on the energy-balance path.
  leaf.wind_speed_ = environment.get_wind_speed();

  // The light chain, from the driver to the leaf. Each stage is named for what
  // it is, because the stages differ by factors of 0.5 and a value handed over
  // at the wrong stage is a plausible number:
  //
  //   SWdown      incoming shortwave, surface, horizontal, above canopy  W m^-2
  //   PPFD_above  = umol_par_per_joule * f_PAR * SWdown                  umol m^-2 s^-1
  //   PPFD_int    = k_I * openness(z) * PPFD_above                       umol m^-2 s^-1
  //               intercepted per unit leaf area: k_I is the extinction
  //               coefficient, which carries the leaf-angle projection
  //   NIR_int     = k_I * openness(z) * (1 - f_PAR) * SWdown             W m^-2
  //
  // The leaf takes the two INTERCEPTED bands and applies its own absorptances
  // (alpha_PAR, alpha_NIR); nothing upstream of it absorbs anything. Both bands
  // see the same openness and k_I, so NIR_int is a fixed multiple of PPFD_int
  // and is formed from it at the leaf rather than integrated over the crown a
  // second time.
  const double f_PAR = environment.f_PAR;
  if (!(f_PAR > 0.0 && f_PAR <= 1.0)) {
    util::stop("TF24: the environment's f_PAR must lie in (0, 1]; got " +
               util::to_string(f_PAR));
  }
  const double PPFD_above =
    phylloptim::umol_par_per_joule * f_PAR * environment.get_SWdown();
  const double NIR_per_PPFD_int =
    (1.0 - f_PAR) / (phylloptim::umol_par_per_joule * f_PAR);

  // Optimise the leaf at a given intercepted PPFD: rebuilds physiology and
  // solves the root-collar water potential, leaving the leaf.* outputs
  // (profit_, transpiration_, soil_consumption_, opt_psi_leaf_, ...) set. Only
  // the light varies between calls; every other input is depth-independent and
  // already computed above.
  //
  // Leaf solves at double, so an active strategy hands it the VALUES of its
  // inputs -- the light, the soil potentials and the conductance are all
  // stripped below -- and record_leaf_outputs puts the rows back on afterwards,
  // by asking the leaf for its own derivatives at the point this solve placed.
  // So the active input has to OUTLIVE the call, which is what `PPFD_int_used`
  // is for: the leaf is handed only its value and cannot give it back. Each leaf
  // output enters the active chain at one place -- the aux stores and
  // leaf.profit_ in net_mass_production_dt, and leaf.soil_consumption_ in
  // evapotranspiration_dt -- so a partial attaches to one expression per
  // output.
  S PPFD_int_used = 0.0;
  // The leaf's inputs at one light, soil and conductance: everything else it
  // reads is the environment's, the same at every call.
  auto seat_leaf = [&](double PPFD, const std::vector<double>& psi,
                       double conductance) -> void {
    leaf.set_physiology(root_network_, PPFD, psi, soil_depths_, conductance,
                        environment.get_atm_vpd(), environment.get_ca(),
                        environment.get_leaf_temp(),
                        environment.get_atm_o2_kpa(),
                        environment.get_atm_kpa());
  };
  auto optimise_at = [&](const S& PPFD_int) -> void {
    PPFD_int_used = PPFD_int;
    // Read on the energy-balance path only; record_leaf_outputs gives it the
    // light's row through the pack's NIR slot.
    leaf.NIR_int_ = NIR_per_PPFD_int * odelia::util::to_passive(PPFD_int);
    // The leaf takes double, and everything the environment supplies already is
    // one. Three of its arguments are this strategy's own and carry a derivative;
    // psi_soil is the one that needs a buffer to be stripped into.
    if constexpr (std::is_same_v<S, double>) {
      seat_leaf(PPFD_int, psi_soil, leaf_specific_conductance_max);
    } else {
      using odelia::util::to_passive;
      psi_soil_value_.resize(psi_soil.size());
      for (size_t a = 0; a < psi_soil.size(); ++a) {
        psi_soil_value_[a] = to_passive(psi_soil[a]);
      }
      seat_leaf(to_passive(PPFD_int), psi_soil_value_,
                to_passive(leaf_specific_conductance_max));
    }
    solve_leaf();
  };

  // Convert canopy openness (0-1) into the intercepted PPFD per unit leaf area:
  // PPFD_above times openness times the extinction coefficient pars.k_I. The
  // light floor (1e-4) matches the mean light's.
  auto PPFD_int_at = [&](const S& light) -> S {
    // Counted, not merely applied. Where this binds the cohort's light is a
    // constant with respect to every other cohort's height, so the row is
    // severed by the guard rather than by the model -- and the field is smooth
    // underneath. An uncounted severance is indistinguishable from a true zero.
    if (light < S(0.0001)) {
      // Below the floor the light this cohort receives stops depending on
      // any other cohort's height. The row is therefore exactly zero for the
      // model AS EVALUATED -- every light below the floor gives a bit-identical
      // census -- so the zero is the derivative of the function on the tape and
      // not an answer withheld. What makes it readable rather than silent is the
      // count, which is taken on this path because it is the only one a clamp
      // severs anything on.
      note_clamp(CLAMP_LIGHT_FLOOR);
    }
    return pars.k_I * std::max(light, S(0.0001)) * PPFD_above;
  };

  // Aggregate the leaf submodel over the crown according to the shading model.
  // The expensive hydraulic optimisation is the unit of work here, so the model
  // choice is about how many times it runs and on what light:
  //  - crown-centre:  one optimisation at the crown-centre light.
  //  - mean-light:    one optimisation at the leaf-area-weighted mean light
  //                   (TF24's established default).
  //  - deep-crown:    one optimisation per crown-depth quadrature point, with
  //                   every leaf output integrated to a leaf-area-weighted mean.
  if (shading_model_ == ShadingModel::CrownCentre) {
    optimise_at(PPFD_int_at(environment.get_environment_at_height(height * eta_c)));
  } else if (shading_model_ == ShadingModel::MeanLight) {
    // Leaf-area-weighted mean canopy openness = integral of (light * q) over the
    // crown (q integrates to one). PPFD_int_at then applies pars.k_I * PPFD_above,
    // exactly reproducing TF24's established average_radiation.
    //
    // Each point's light is floored at 1e-4 first (the original rationale was
    // never recorded; preserved as-is). Counted here, because this is where the
    // floor binds: the crown weights sum to one, so a mean of floored values
    // cannot fall below the floor itself.
    const auto floored = [this](const auto& light) {
      using X = std::decay_t<decltype(light)>;
      if (light < X(0.0001)) {
        note_clamp(CLAMP_LIGHT_FLOOR_CROWN);
        return X(0.0001);
      }
      return X(light);
    };
    optimise_at(PPFD_int_at(
        environment.light_availability.crown_integral(height, crown, floored)));
  } else { // DeepCrown
    if constexpr (std::is_same_v<S, double>) {
      const std::vector<double> nodes =
        function_integrator.integrate_vector_x(0.0, height);
      const size_t nn = nodes.size();
      std::vector<double> profit_y(nn), trans_y(nn), eup_y(nn), psi_y(nn),
        root_psi_y(nn), gco2_y(nn), assim_y(nn), tleaf_y(nn);
      std::vector<std::vector<double>> soil_y(
        soil_number_of_depths_, std::vector<double>(nn));
      for (size_t i = 0; i < nn; ++i) {
        const S qi = canopy_shape.q_from_height(nodes[i], height);
        optimise_at(PPFD_int_at(environment.get_environment_at_height(nodes[i])));
        profit_y[i]   = leaf.profit_ * qi;
        trans_y[i]    = leaf.transpiration_ * qi;
        eup_y[i]      = leaf.E_up_ * qi;
        psi_y[i]      = leaf.opt_psi_leaf_ * qi;
        root_psi_y[i] = leaf.opt_root_psi_ * qi;
        gco2_y[i]     = leaf.stom_cond_CO2_ * qi;
        assim_y[i]    = leaf.assim_colimited_ * qi;
        // Leaf temperature varies through the crown because the light does, so it
        // must be integrated like every other leaf output. Left out, `Tleaf_`
        // would report whichever node the loop happened to end on -- and that node
        // is neither the crown top nor the centre.
        tleaf_y[i]    = leaf.Tleaf_ * qi;
        for (int a = 0; a < soil_number_of_depths_; ++a) {
          soil_y[a][i] = leaf.soil_consumption_[a] * qi;
        }
      }
      // Integrate each leaf output to its leaf-area-weighted crown mean (q
      // integrates to one over the crown). soil_consumption_ feeds the patch
      // water balance, so it must be the depth-integrated total; the rest are
      // diagnostics reported through compute_rates.
      leaf.profit_          = function_integrator.integrate_vector(profit_y, 0.0, height);
      leaf.transpiration_   = function_integrator.integrate_vector(trans_y, 0.0, height);
      leaf.E_up_            = function_integrator.integrate_vector(eup_y, 0.0, height);
      leaf.opt_psi_leaf_    = function_integrator.integrate_vector(psi_y, 0.0, height);
      leaf.opt_root_psi_    = function_integrator.integrate_vector(root_psi_y, 0.0, height);
      leaf.stom_cond_CO2_   = function_integrator.integrate_vector(gco2_y, 0.0, height);
      leaf.assim_colimited_ = function_integrator.integrate_vector(assim_y, 0.0, height);
      leaf.Tleaf_           = function_integrator.integrate_vector(tleaf_y, 0.0, height);
      for (int a = 0; a < soil_number_of_depths_; ++a) {
        leaf.soil_consumption_[a] =
          function_integrator.integrate_vector(soil_y[a], 0.0, height);
      }
    } else {
      util::stop("shading_model '" + this->control.shading_model +
                 "' is not differentiable: its crown means pass through Leaf, "
                 "which carries double. Use crown-centre or mean-light "
                 "shading for a gradient.");
    }
  }


  //TODO: one point constant ratio and integral width for daylength
  // convert assimilation per leaf area per second (umol m^-2 s^-1) to canopy-level total yearly assimilation (mol yr^-1)
  // converts to canopy area, then years, then mols
  //
  // ⚠️ GROWTH IS BILLED ON THE CARBON KEPT, phylloptim's `carbon()` (the
  // objective plus `shadow_cost`), NOT ON THE OBJECTIVE. On TF24_floor the objective deducts `lambda_o * E`, the shadow
  // price of water: the value of the water in its best alternative use, which for
  // a leaf is assimilation later. No carbon is lost when the plant pays it; it
  // changes the aperture chosen and nothing else, which is what a Lagrange
  // multiplier does. Feeding the objective straight into growth would tax the
  // plant by carbon it never spent: measured on a 5 m plant at PPFD 1800 and theta
  // 0.25, the objective understates the carbon kept by 2.3% at lambda_o = 1e4,
  // 9.4% at 5e4, 15.4% at 1e5 and 22.6% at 2e5.
  //
  // On TF24_floor `carbon()` is `profit_ + lambda_o * transpiration_`, bit for
  // bit what this summed before phylloptim billed it, and `profit_` itself at the
  // default price. No second canopy integral is needed: it reads the stored
  // `profit_` and `transpiration_`, which the DeepCrown branch integrated against
  // the same weights, and it is linear in both.
  S profit_ = leaf.carbon();
  if constexpr (!std::is_same_v<S, double>) {
    record_leaf_outputs(PPFD_int_used, psi_soil,
                        leaf_specific_conductance_max, NIR_per_PPFD_int);
    profit_ = leaf_carbon_;
  }

  // --- How much would more conducting area be worth? (#516) ------------------
  // d(profit)/d(kmax) at the operating point the leaf just found.
  //
  // ⚠️ THIS COSTS TWO FURTHER LEAF SOLVES, and there is no way around it. The
  // obvious shortcut is the envelope theorem: profit is already optimised, so
  // re-evaluate at the collar potential already found and the indirect term
  // vanishes. THAT ARGUMENT IS WRONG HERE. It holds for opt_psi_leaf_, which
  // really is chosen to maximise profit. It does NOT hold for opt_root_psi_,
  // which is found by a ROOT FIND -- the collar potential at which the soil-root
  // network's supply matches the leaf's demand. A constraint, not an optimum, so
  // d(profit)/d(collar) is not zero and dropping the indirect term drops a real
  // quantity. Measured: holding the collar fixed under-reports d(profit)/d(kmax)
  // by 15-19 per cent, which put the controller's zero at 1.22x the pipe-model
  // ratio where growth actually peaks at 1.28x. Re-solving matches a finite
  // difference of the model's own assimilation to 0.09 per cent.
  //
  // THE STEP IS MEASURED, NOT ASSUMED. profit comes through a solve with its own
  // internal tolerances, so its effective precision is far short of machine
  // epsilon. 1e-3 relative sits clear of that noise floor while staying deep in
  // the linear regime; sweeping it 1e-6 to 1e-3 moved the result by 0.3 per
  // cent, so the step is not what limits accuracy here.
  //
  // Placed after record_leaf_outputs, which reads the leaf at the operating point
  // the profit came from; the third solve then puts it back there for the aux
  // reads in compute_rates. The difference is taken in double and carries no
  // rows, so a gradient through it would be missing the controller's own
  // sensitivity: refused rather than returned. Skipped entirely at a_sw = 0, the
  // default, so the hot path is untouched unless sapwood acclimation is on.
  // (prepare_strategy refuses a_sw > 0 under deep-crown, whose profit is a crown
  // integral rather than one radiation's.)
  dprofit_dkmax_ = 0.0;
  if (odelia::util::to_passive(pars.a_sw) > 0.0) {
    using odelia::util::to_passive;
    if constexpr (!std::is_same_v<S, double>) {
      refuse("sapwood acclimation (a_sw > 0) is not differentiable: its "
             "sensitivity is a finite difference through the leaf, which "
             "carries no rows");
    }
    const double kmax = to_passive(leaf_specific_conductance_max);
    const double profit_base = leaf.carbon();
    const double dk = 1e-3 * kmax;
    psi_soil_value_.resize(psi_soil.size());
    for (size_t a = 0; a < psi_soil.size(); ++a) {
      psi_soil_value_[a] = to_passive(psi_soil[a]);
    }
    auto solve_at = [&](double conductance) {
      seat_leaf(to_passive(PPFD_int_used), psi_soil_value_, conductance);
      solve_leaf();
    };
    solve_at(kmax + dk);
    dprofit_dkmax_ = (leaf.carbon() - profit_base) / dk;
    // Put the leaf back where it was. Every aux compute_rates writes reads the
    // leaf's members, and a leaf left at the perturbed conductance would report
    // a plant that does not exist.
    solve_at(kmax);
  }

  const S assimilation_ = profit_ * area_leaf_* 60*60*12*365/1e6;
  // const double assimilation_ = assimilation(environment, height, area_leaf_);
  const S mass_coarse_root_ = mass_coarse_root(mass_sapwood_);
  const S respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_,
                mass_coarse_root_);

  // --- Is the marginal leaf paying for itself? (#516) ------------------------
  // The gate reads this rather than the reserve pool, following Manzoni et al.
  // 2015: maximising a season's carbon gain over leaf area gives dG/dL = 0, i.e.
  // shed once a leaf stops covering its own upkeep. That criterion fires while
  // the plant is still carbon-rich -- waiting for reserves to run down IS the
  // evergreen strategy in that paper, and it is the one that loses a long
  // drought. Robinson's mulga agree from the field: 85 per cent of trees shed,
  // the healthiest included, and shedders are indistinguishable from survivors
  // until a year in.
  //
  // ⚠️ ONLY COSTS THAT SCALE WITH LEAF AREA BELONG HERE, because only those are
  // what a marginal leaf brings with it: its own upkeep, its fine root and its
  // bark (both pinned to leaf area). Sapwood is excluded deliberately -- it is
  // carried whether or not this leaf exists, so charging it here would make the
  // marginal leaf look unprofitable for a cost it is not responsible for, and
  // would make a tall plant shed for being tall.
  const S scaling_respiration_ =
    respiration_leaf(mass_leaf_) + respiration_bark(mass_bark_) +
    respiration_root(mass_root_);
  const S scaling_turnover_ =
    turnover_leaf(mass_leaf_) + turnover_bark(mass_bark_) +
    turnover_root(mass_root_);
  const S scaling_cost_ =
    pars.a_bio * pars.a_y * scaling_respiration_ + scaling_turnover_;
  const S scaling_gain_ = pars.a_bio * pars.a_y * assimilation_;
  // Dimensionless, and zero exactly at Manzoni's break-even, so the gate's
  // centre is a pure number whose default of 0 IS that optimum.
  marginal_leaf_return_ = scaling_cost_ > 0.0
    ? S((scaling_gain_ - scaling_cost_) / scaling_cost_)
    : S(0.0);

  // The gate reads the TRACKED marginal return, not the instantaneous one just
  // computed: assimilation swings with the weather, and a plant that shed on
  // every dry week would be neither physical nor integrable.
  const S replacement_leaf =
    1.0 - (1.0 - replacement_fraction(tracked_marginal)) *
            shedding_gate(leaf_departure);
  const S replacement_sapwood =
    sapwood_replacement_fraction(replacement_leaf, sapwood_departure);
  withheld_leaf_    = 1.0 - replacement_leaf;
  withheld_sapwood_ = 1.0 - replacement_sapwood;

  const S turnover_ =
    turnover(replacement_leaf, replacement_sapwood, mass_leaf_, mass_bark_,
             mass_sapwood_, mass_root_, mass_coarse_root_);
  // The plant's maintenance bill, which the sapwood controller divides by to
  // get a dimensionless return. Used instead of assimilation or profit because
  // it is STRICTLY POSITIVE while the plant has tissue: profit goes negative in
  // exactly the drought where the controller most needs a sign, and a
  // denominator that changes sign would flip the direction the stem moves.
  maintenance_flux_ = pars.a_bio * pars.a_y * respiration_ + turnover_;
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// Base TF24: place the operating point the run found at this rate evaluation, or
// optimise the root-collar water potential where it kept none -- a branch that exited
// on feasibility, or an evaluation the run addressed no record against.
//
// ⚠️ THE TRY/CATCH TURNS AN UNPHYSICAL PSI PROBE INTO A REJECTED STEP rather than
// a dead run. phylloptim raises its own infeasible_error, which derives from
// std::runtime_error and NOT from odelia::util::DomainError -- they are siblings
// -- so odelia's stepper cannot recognise it, and the throw kills the whole solve
// having taken zero steps (#608 measured exactly this). Translating it lets odelia
// shrink and retry, and it stops only if the minimum step still cannot reach a
// feasible probe, reporting phylloptim's own message when it does.
//
// Deliberately narrow: only infeasible_error is translated. A util::stop() from
// phylloptim, or any other exception, still propagates, so a bug stays a bug
// instead of becoming step-shrinking until "Cannot achieve the desired accuracy".
template <typename S>
void TF24_Strategy<S>::solve_leaf() {
  const leaf_solved_point recorded = leaf_points->load();
  try {
    if (recorded.kind == Leaf::OperatingPointKind::Unsolved) {
      // THE CURVE-TYPED FORM, which names the curve the surface below reads.
      // (find_root_collar_psi() now solves the seated curve too.)
      leaf.template find_root_collar_psi_for<
          phylloptim::Leaf::CostCurve::TF24_floor>();
    } else {
      // ⚠️ ONE CALL, because the collar and the arm have to go back together:
      // evaluating at a target restores every number and then tags the point
      // `prescribed`, and operating_collar switches on the kind.
      leaf.replay_operating_point(recorded.collar, recorded.kind);
    }
  } catch (const phylloptim::util::infeasible_error& e) {
    odelia::util::stop_domain(std::string("leaf solve infeasible: ") + e.what());
  }
  leaf_points->store({leaf.opt_root_psi_, leaf.operating_point_kind()});
  // The classification is decided by the branch taken and then overwritten by
  // the next plant, so without a tally the only route to its incidence is a
  // refusal message -- which reports the FIRST non-interior point and nothing
  // about how many followed it or what kinds they were.
  ++operating_point_counts[
      static_cast<size_t>(leaf.operating_point_kind())];
}

}

#endif
