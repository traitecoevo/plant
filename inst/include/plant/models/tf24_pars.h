// -*-c++-*-
#ifndef PLANT_PLANT_TF24_PARS_H_
#define PLANT_PLANT_TF24_PARS_H_

#include <plant/pars_table.h>
#include <array>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <vector>

namespace plant {

// Biological (user-settable) parameters for the TF24 model. Held as a value
// member `pars` on TF24_Model and exposed to R as a nested RcppR6 list class
// (access as `s$pars$lma`). Only R-settable parameters live here; derived
// quantities (eta_c, height_0, ...), the embedded Leaf model, solver
// tolerances and hard-coded hydraulic-root constants stay as plain members
// on the strategy.
template <typename S = double>
struct TF24_Pars {
  using value_type = S;

  // A default member initialiser has no block to declare the library pow in, so
  // the derived defaults below raise their base through here instead.
  static S power(const S& base, const S& exponent) {
    using std::pow;
    return pow(base, exponent);
  }

  // * Core traits
  S lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  S rho       = 608.0;      // Wood density [kg/m3]
  S hmat      = 16.5958691; // Height at maturation [m]
  S omega     = 3.8e-5;     // Seed mass [kg]
  // * Individual allometry
  S eta       = 12.0;       // Canopy shape parameter [dimensionless]
  S theta     = 1.0/4669;   // Sapwood area per leaf area [dimensionless]
  S a_l1      = 5.44;       // height with 1m2 leaf [m]
  S a_l2      = 0.306;      // scaling of height with leaf area
  S a_r1      = 0.07;       // Fine-root mass per leaf area [kg / m]
  // Coarse (structural) root mass per unit sapwood mass. Coarse roots are
  // treated as the below-ground continuation of the sapwood cylinder, so they
  // inherit its height scaling and the structural root:shoot ratio is
  // size-invariant. Zero by default, so the pool is absent unless switched on.
  // Only a_r1 feeds the depth-distributed water uptake -- that is a fine-root
  // process, and coarse roots must stay out of it (see #349).
  S a_cr1     = 0.0;        // [dimensionless]
  S a_b1      = 0.17;       // Ratio of bark area : sapwood area
  // * Production
  S r_s    = 4012.0 / 608.0; // Sapwood respiration per stem mass
  S r_b    = 2.0 * r_s;      // Bark respiration (assumed 2 x sapwood)
  S r_r    = 217.0;          // Fine-root respiration per mass
  // Coarse roots are woody, so this defaults to the sapwood rate. As with r_b,
  // a default and not a link: setting r_s afterwards does not move it.
  S r_cr   = r_s;            // Coarse-root respiration per mass
  S r_l    = 39.27 / 0.1978791; // Leaf dark respiration per leaf mass
  S a_y    = 0.7;            // Carbon conversion parameter
  S a_bio  = 2.45e-2;        // CO2 -> dry mass [kg / mol]
  S k_l    = 0.4565855;      // Leaf turnover [/yr]
  S k_b    = 0.2;            // Bark turnover [/yr]
  S k_s    = 0.2;            // Sapwood turnover [/yr]
  S k_r    = 1.0;            // Fine-root turnover [/yr]
  S k_cr   = k_s;            // Coarse-root turnover [/yr] (woody)
  S a_p1   = 151.177775377968;   // LRC hyperbola [mol CO2 / yr / m2]
  S a_p2   = 0.204716166503633;  // LRC hyperbola shape
  // * Seed production
  S a_f3   = 3.0 *  3.8e-5;  // Accessory cost of reproduction [kg/seed]
  S a_f1   = 1.0;            // Maximum allocation to reproduction
  S a_f2   = 50;             // Size range across which individuals mature
  // * Mortality parameters
  S S_D    = 0.25;           // Probability of survival during dispersal
  S a_d0   = 0.1;            // Parameter for seedling survival
  S d_I    = 0.01;           // Baseline intrinsic mortality [/yr]
  // a_dG1 / a_dG2 shape the *storage-dependent* growth mortality:
  //   mortality_storage_dependent_dt(r) = a_dG1 * exp(-a_dG2 * r),  r = S/S_max.
  // ⚠️ DO NOT feed productivity_area to this exponent in place of r: under deep
  // carbon deficit it overflows to ~1e32. Reserves keep r in [0,1], so the term
  // is bounded in [a_dG1*e^-a_dG2, a_dG1].
  S a_dG1  = 5.5;            // Max growth-related (low-reserve) mortality [/yr]
  S a_dG2  = 20.0;           // Sensitivity of mortality to relative reserves
  // * NSC storage pool (#517) -- buffers growth & mortality against short-term
  //   productivity swings. See Stefaniak et al. 2026 (plantNSC) for the design.
  //   Carbon first charges storage; growth/reproduction are mobilised from it,
  //   but gated so growth only proceeds when reserves are ample (a plant should
  //   not grow down its stores). Reserves-based mortality rises as they deplete.
  S a_st1  = 0.10;           // Storage capacity per unit sapwood mass [kg NSC / kg]
  S a_st2  = 0.1;            // Reserve fraction at which growth is half-on [0-1]
  S a_st3  = 0.8;            // Initial storage at birth [fraction of capacity]
  // Added to the storage pool's own relaxation time, storage_max/(charge +
  // drain) [yr].
  S storage_relaxation_offset = 7.0 / 365.0;
  // * Flexible allometry (#516) -- a plant under carbon stress declines to
  //   replace the leaf and sapwood its turnover destroys, so its canopy thins
  //   instead of its mortality spiking, and rebuilds when carbon returns.
  //   Settable as traits, being TF24_Pars fields; not hyperparameterised,
  //   there being no trade-off to encode until #512 reports one.
  // a_pl0 = 0 is the fixed allometry EXACTLY, not approximately: it zeroes the
  // withheld fraction, which leaves both departures resting at zero.
  S a_pl0  = 0.0;            // Max fraction of leaf/sapwood turnover left unreplaced [0-1]
  // ⚠️ a_pl1 GATES ON THE MARGINAL LEAF'S CARBON BALANCE, NOT ON RESERVES. It is
  // the dimensionless return (gain - cost)/cost at which withholding is half its
  // maximum, so 0 is exactly Manzoni et al. (2015)'s optimum -- shed once a leaf
  // stops covering its own upkeep. An earlier version gated on the NSC reserve
  // fraction; that fires only after the pool has drained, which is the evergreen
  // strategy Manzoni shows losing a long drought, and it was measured here to buy
  // a factor of 1.0-1.13 on survivorship. Negative values make a plant hold its
  // canopy past break-even; positive values make it shed while still in profit.
  S a_pl1  = 0.0;            // Marginal leaf return at which withholding is half its max
  S a_pl2  = 0.2;            // Log-gap scale for rebuilding, and for the canopy floor
  // Smallest canopy a plant will thin to, as a fraction of what its height
  // implies. Two jobs. Biologically it is the evergreen habit: Robinson's mulga
  // that de-leafed severely died, so for an evergreen the floor is real rather
  // than a numerical convenience. Numerically it is what keeps the state in a
  // range where the arithmetic means something -- without it the canopy decays
  // exponentially with no bound, reaching 2e-8 for a fast-leaved species and
  // 0.000 for a suppressed cohort in an ordinary wet stand, with divisions by
  // leaf area reaching 1e9 inside a solver step shared with live cohorts.
  //
  // Set it very small for a drought-deciduous species, which sheds to nothing and
  // re-flushes. That is supported and should stay reachable, but it is not the
  // habit this model is aimed at, and the departure coordinate degrades as the
  // floor approaches zero (phi -> -inf).
  S a_pl3  = 0.05;           // Minimum canopy, as a fraction of the preferred [0-1]
  // Rate at which the tracked marginal return chases the instantaneous one, /yr.
  // The gate must read a RUNNING MEAN: assimilation swings with the weather, and
  // a plant that shed and re-flushed on every dry week would be neither
  // physical nor integrable. 1/a_pl4 is the averaging window, so the default is
  // about a four-month memory -- long against weather, short against a drought.
  S a_pl4  = 3.0;            // Acclimation rate of the tracked marginal return [/yr]
  // Rate at which sapwood per leaf area chases the ratio that maximises GROWTH,
  // /yr. Zero (the default) pins it to the pipe model and leaves every result
  // unchanged.
  //
  // The controller integrates the marginal return, and an integrator IS an
  // averager: a small rate makes the stem respond to the mean of a signal that
  // swings with the weather, without needing a second state to track it.
  //
  // ⚠️ The objective is growth, NOT net production. Priced properly, extra stem
  // is paid for in forgone leaf-area growth, and the two optima differ a lot:
  // production peaks near 3.3x the pipe-model ratio, growth near 1.4-1.8x.
  S a_sw   = 0.0;            // Sapwood acclimation rate [/yr]
  // * Light capture
  // The canopy extinction coefficient, k = G / cos(theta): 0.5 is a spherical
  // leaf-angle distribution under an overhead sun. It does TWO jobs. It is the
  // Beer's-law exponent in the shading this plant casts on others
  // (compute_competition), and it is the factor that turns irradiance on a
  // horizontal plane into irradiance INTERCEPTED per unit leaf area -- the
  // leaf-angle projection -- in the light handed to the leaf. Applied once
  // there; the leaf applies the absorptance.
  S k_I = 0.5;
  // * Leaf hydraulic / photosynthesis traits (default Eucalyptus saligna)
  S vcmax_25 = 96;
  // ⚠️ THE LEAF TRAITS CARRY PHYLLOPTIM'S NAMES, and that is a rule rather than a
  // tidy-up. These are handed to the `Leaf` constructor POSITIONALLY, so a
  // mismatch between what a name says here and which slot it lands in is silent:
  // `g1_TF24` was the old name for what phylloptim calls `TF24_cost_scale`, and it
  // reads as a stomatal-slope g1 -- a quantity the leaf also has, under `g1`, in
  // the Medlyn path. One name per quantity, spelled the same on both sides of the
  // boundary, is what makes the hand-over at prepare_strategy() checkable by eye.
  //
  // ⚠️ `stem_b` AND `psi_crit` ARE DERIVED HERE AND HANDED OVER TO NOBODY.
  // phylloptim takes (stem_P50, stem_c) and derives both itself, by these same
  // formulas -- so these two members are the trade-off function's reporting
  // copies, not inputs to the model. Setting one does not change the curve.
  // stem_P50 and stem_c are the frozen TF24 trade-offs' values at the default
  // K_s (frozen_TF24_hyperpar in R/default_parameters.R), so a
  // strategy built from traits and the bare default share one curve (#611).
  S stem_P50 = 2.8887256653360196;
  // Sapwood-specific conductivity of the TERMINAL segment: the whole-stem value
  // 1 of the height-linear model, back-derived to the tip by
  // TF24_K_s_from_whole_stem(1), so that the two models' resistances agree at
  // TF24_H_ANCHOR = 1 m. The ratio is 2.9782.
  //
  // Only theta/K_s is identifiable, so agreement with the height-linear model
  // holds at exactly ONE height and nowhere else. That is by construction:
  // changing the height dependence is the object of the exercise.
  S K_s = 0.33577377016801868;
  S stem_c = 2.04;
  S stem_b = stem_P50 / power(-log(1 - 50.0 / 100.0), 1 / stem_c);
  S psi_crit = stem_b*power(log(1/0.05),1/stem_c); // derived from stem_b and stem_c
  // Read by nothing: not passed to the leaf, not read by any equation here. Kept
  // only because it is in the R-visible parameter list.
  S beta1 = 20000;
  S TF24_beta2 = 1.5;
  S TF24_cost_scale = 7.5;
  // THE PRICE OF WATER AS TRANSPIRATION GOES TO ZERO, umol C (kg H2O)^-1: the
  // `lambda_o` of phylloptim's TF24_floor cost curve, which TF24 and TF24f are
  // seated on (#634).
  //
  // WHY IT EXISTS. Every conductance-loss cost -- TF24's included -- prices water
  // at zero as E goes to zero, since no conductivity is lost when nothing flows,
  // so water is free precisely when it is abundant. TF24_floor splits the cost
  // into a part depending on potential alone and a part linear in the flux,
  // Theta(E) = Theta~(psi) + lambda_o*E, with Theta~ being TF24's own cost at
  // TF24's own traits. So this is the curve's ONLY parameter, and **TF24 is
  // TF24_floor at lambda_o = 0, at identical parameter values** -- which makes
  // "is TF24 missing a price of water?" a one-restriction question.
  //
  // ⚠️ DEFAULTS TO ZERO, WHICH IS TF24 EXACTLY -- not approximately. phylloptim
  // asserts the reduction bit-for-bit rather than to a tolerance, because each
  // term is the parent's own expression so zeroing the price adds an exact zero to
  // TF24's exact value (it survives fused multiply-add for the same reason). So
  // every result this repo has recorded is unchanged at the default, and a moved
  // regression baseline means a mistake rather than a new expected value.
  //
  // ⚠️ IT IS A SHADOW PRICE, NOT A COST, and that distinction is why
  // net_mass_production_dt books phylloptim's `carbon()`, the objective with the
  // price added back, before converting to carbon. `Theta~` is carbon actually forgone; `lambda_o` is the
  // value of water in its best alternative use, which for a leaf is assimilation
  // later. Paying it changes the aperture chosen and loses no carbon.
  //
  // ⚠️ ITS SCALE IS SET BY THE LEAF, not chosen freely: phylloptim's own marginal
  // cost of water runs ~9e4 to 3e5 in these units at its defaults, so a value far
  // outside that band pins the optimum against a bracket bound and the answer
  // describes the bracket rather than the model.
  S TF24_floor_lambda_o = 0.0;
  S jmax_25 = vcmax_25*1.64;
  // Quantum yield of electron transport per ABSORBED photon. The leaf absorbs
  // alpha_PAR of the intercepted PPFD, so 0.30 / 0.85 is the old incident-basis
  // yield of 0.30 restated (phylloptim 0.13.0 renamed `a` to `phi_J`).
  S phi_J = 0.30 / 0.85;
  S curv_fact_elec_trans = 0.7;
  S curv_fact_colim = 0.99;
  // Dark respiration at 25 C, the term net assimilation subtracts. It is handed
  // over now that set_traits takes it, rather than left to the leaf's own
  // default -- so a parameter can move it and a gradient can see it.
  S R_d_25 = 1.44;
  S var_sapwood_volume_cost = 1;
  // nitrogen allocation traits (parameterised from Austraits 4.1.0)
  S nmass_l = 13e-3; // kg N kg^-1 mass
  S nmass_s = 1.98e-3; // kg N kg^-1 mass
  S nmass_b = 3.40e-3; // kg N kg^-1 mass
  S nmass_r = 3.35e-3; // kg N kg^-1 mass
  S dmass_dN = 0; // change in mass per change in kg kg^-1 N
  // shape exponent for the Q() root-fraction-with-depth profile
  S root_depth_shape_eta = 0.2;
  // * Root hydraulics
  // Root vulnerability curve, proportion of conductivity =
  // exp(-(psi/root_b)^root_c). Declared before root_psi_crit, which derives
  // from them. Both are settable from R, because the defaults put root shutoff
  // at ~5.87 MPa: too conservative for taxa that operate below that (e.g.
  // Acacia aneura).
  S root_c = 2.680147;
  // THE SCALE PARAMETER IS P50, matching the stem side and phylloptim, which
  // takes (root_c, root_P50) and derives the other two itself by these same
  // formulas. 3.4 reproduces the old root_b = 3.898245 to ten significant
  // figures: derived it is 3.8982451221145307, +3.1e-10 relative, and
  // root_psi_crit moves with it by the same fraction.
  S root_P50 = 3.4;
  S root_b = root_P50 / power(-log(1 - 50.0 / 100.0), 1 / root_c);
  // Potential at 5% remaining root conductivity [MPa]. Derived, exactly as
  // psi_crit is from b and c: if you set root_b or root_c directly, set this
  // too, or the vulnerability curve and the shutoff threshold disagree.
  S root_psi_crit = root_b*power(log(1/0.05),1/root_c);
  // Maximum rooting depth [m]. Rooting depth is min(height, rooting_depth_max),
  // so this also bounds the depth over which roots can draw water. Should not
  // exceed the soil column depth (TF24_Environment `depth`, default 1.5 m):
  // layers below the column do not exist, so deepening roots alone gains
  // nothing without deepening the soil as well.
  S rooting_depth_max = 1.5;
  // * Stem hydraulic path
  // Within-plant anatomical profiles along the flow path, parameterised by
  // distance from the apex L and anchored at the terminal segment. See
  // plant/stem_hydraulics.h for the closed form and what it assumes.
  //
  // Setting all three to zero collapses the path integral to a resistance
  // linear in height -- the model before this existed -- bit for bit, which is
  // what the regression test asserts. theta_c is zero by default, so `theta`
  // keeps its whole-plant meaning and no parameter-file migration is needed
  // yet.
  //
  // Conduit widening exponent: D(L) = D_tip*(L/L_tip)^D_c. Reaches the model
  // only through beta = 2*D_c + theta_c -- the diameter itself is never
  // evaluated, and D_tip is not a parameter until the sec. 4.2 diagnostic.
  // Named D_c, not b, because b is already the Weibull vulnerability scale
  // above; `_c` means "exponent" here, as in c and root_c.
  // Measured, not fitted: conserved across terrestrial vascular plants, with a
  // within-stem range of roughly 0.1-0.3.
  S D_c = 0.2;
  // Huber-profile exponent: theta(L) = theta*(L/L_tip)^(-theta_c). NOTE THE
  // MINUS SIGN. theta falls basipetally while the Huber value 1/theta rises, so
  // a positive theta_c means less leaf area supported per unit sapwood towards
  // the base -- that is the compensation mechanism.
  //
  // NOT YET IMPLEMENTED: prepare_strategy() throws on any non-zero value. theta
  // is not a hydraulics-only trait -- it also sets area_sapwood, area_bark,
  // mass_sapwood (hence construction cost, respiration, turnover and NSC
  // capacity) and the hard-coded dmass_sapwood_darea_leaf derivative, all of
  // which still read a flat pars.theta. Profiling it on the hydraulic side
  // alone would give a plant that conducts as though theta varied and is built
  // as though it did not. The two uses are the same trait and must move
  // together, so the field is declared and refused rather than half-applied.
  //
  // Name clash to be aware of: phylloptim's Leaf spells soil water content
  // theta_, theta_w_ and theta_fc_ (m^3 m^-3), all R-visible, so s$pars$theta_c
  // (dimensionless) and leaf$theta_ coexist in one session.
  S theta_c = 0.0;
  // Terminal segment length [m]. Not an innocuous numerical cutoff: it enters
  // the resistance with elasticity beta ~ 0.6 and trades off exactly against
  // both K_s and theta, which are identifiable only in the grouping
  // theta*L_tip^beta/K_s. It must therefore come from the SAME terminal-segment
  // definition over which K_s and theta were measured -- none of the three may
  // be calibrated independently of the others (invariance criterion I5).
  S L_tip = 0.02;
  // Germination
  S recruitment_decay = 0.0;
  // Penman-Monteith leaf energy balance (#523). use_energy_balance gates PM
  // (0 = off, today's Tleaf=Tair behaviour; != 0 = on); default off preserves
  // backward compatibility. d is the characteristic leaf dimension (m) for the
  // aerodynamic resistance ra = C_ra*sqrt(d/U0); inert while PM is off.
  S use_energy_balance = 0.0;
  S d = 0.05;
  // Leaf absorptances of intercepted PAR and near-infrared, handed to the leaf.
  // alpha_PAR reaches photosynthesis (with phi_J) and the energy balance;
  // alpha_NIR only the energy balance, so it is inert while PM is off.
  S alpha_PAR = 0.85;
  S alpha_NIR = 0.20;
  // Crown leaf area index: the leaf area a crown carries per unit of ground it
  // covers, so crown area = area_leaf / L_crown. Read by TF24w alone (#681),
  // where it sets the footprint of a plant's private root-zone soil column;
  // TF24 reads nothing from it. Default to be set from crown allometry; 2 is a
  // placeholder.
  S L_crown = 2.0;
  // Residual (minimum) leaf conductance to water at 25 C, mol H2O m^-2 s^-1:
  // phylloptim's `g_min` (phylloptim #182), the flux a leaf loses with its
  // stomata shut. Zero by default, which is exactly the model without it: the
  // leaf tests `g_min == 0.0` and returns exact zeros for the residual terms, so
  // no result moves at the default. Carried here so a parameter set can reach
  // every leaf trait through `pars`; a non-zero value is untested in plant.
  S g_min = 0.0;

  // Where a parameter lives and what it is called. Nothing here says what its
  // gradient will be: the sweep computes one for every parameter except the
  // `undifferentiable` names below.
  using ad_parameter = pars_field<TF24_Pars, S>;
  using no_gradient = plant::no_gradient;

  // The parameters with no gradient column, and why each has none. This is the
  // whole of what the gradient needs declared: everything else in the table gets
  // a column, so a parameter added to the model is differentiable because nobody
  // said it was not. Asking for one of these refuses in the words beside it.
  //
  // ⚠️ A PARAMETER NO EQUATION READS BELONGS HERE, and giving it a column instead
  // is this design's worst failure mode: the sweep returns an exact zero, which
  // is the correct derivative and reads as an answer. It also disarms the guard
  // -- "every registered TF24 parameter reaches an output" -- that catches a
  // parameter added to the model and never wired into it. `a_f3` is the contrast
  // and is NOT here: it reaches fecundity_dt, so it is registered, and its census
  // column is zero only because no census metric reads that rate.
  static constexpr std::array undifferentiable{
      // No row can be recorded.
      no_gradient{"eta",
                  "the canopy profile's u^eta reaches u = 0, where the recorded "
                  "derivative u^eta*log(u) is 0*(-inf) and the guard returns a "
                  "constant zero instead, so a row here would be a silently "
                  "wrong zero rather than a NaN"},
      no_gradient{"root_depth_shape_eta",
                  "Q()'s u^eta_x reaches u = 0 the same way, and the guard "
                  "supplies the cumulative fraction there outright"},
      no_gradient{"d", "the leaf supplies no row for it"},
      no_gradient{"alpha_PAR",
                  "a leaf property the leaf holds as a constant, so it supplies "
                  "no row; phi_J, which it multiplies, carries the column"},
      no_gradient{"alpha_NIR", "the leaf supplies no row for it"},
      no_gradient{"use_energy_balance",
                  "a gate, compared rather than differentiated"},
      // The coarse-root pool is opt-in: at the default a_cr1 = 0 it is empty, so
      // r_cr and k_cr multiply zero and a column would be an exact zero that
      // reads as an answer. No gradient fixture runs with the pool on.
      no_gradient{"a_cr1",
                  "the coarse-root pool is opt-in and has no gradient fixture"},
      no_gradient{"r_cr", "the coarse-root pool is empty at the default a_cr1 = 0"},
      no_gradient{"k_cr", "the coarse-root pool is empty at the default a_cr1 = 0"},
      // Flexible allometry (#516) is opt-in the same way: at a_pl0 = 0 nothing
      // is withheld, so the other four gate parameters multiply zero.
      no_gradient{"a_pl0",
                  "flexible allometry is opt-in and has no gradient fixture"},
      no_gradient{"a_pl1", "the replacement gate is shut at the default a_pl0 = 0"},
      no_gradient{"a_pl2", "the replacement gate is shut at the default a_pl0 = 0"},
      no_gradient{"a_pl3", "the replacement gate is shut at the default a_pl0 = 0"},
      no_gradient{"a_pl4", "the replacement gate is shut at the default a_pl0 = 0"},
      no_gradient{"a_sw",
                  "sapwood acclimation's sensitivity is a finite difference "
                  "through two further leaf solves, which carry no rows, so a "
                  "run with it on refuses the gradient; at the default a_sw = 0 "
                  "its column would be an exact zero"},
      no_gradient{"L_crown",
                  "read by TF24w alone, which refuses the gradient path: its "
                  "private root-zone water couples cohorts through the patch, "
                  "outside the per-cohort block the recording describes"},
      // The residual conductance is opt-in the same way: at the default
      // g_min = 0 the leaf's residual terms are exact zeros by a test on the
      // value, so a column would be an exact zero that reads as an answer.
      no_gradient{"g_min",
                  "the residual conductance is opt-in and has no gradient "
                  "fixture; at the default g_min = 0 its column would be an "
                  "exact zero"},
      // Refused by the model, so a row would describe half of it.
      no_gradient{"theta_c",
                  "prepare_strategy refuses any non-zero value, because theta "
                  "also sets sapwood and bark area, construction cost, "
                  "respiration and storage capacity and those still read a flat "
                  "theta -- the path length does move with it, so a column here "
                  "would be a live number for the hydraulic half of a trait the "
                  "carbon budget does not yet follow"},
      // No equation reads them.
      no_gradient{"stem_b", "derived from (stem_P50, stem_c); phylloptim derives it itself, so setting it reaches nothing and a row here would describe a path the model does not take"},
      no_gradient{"psi_crit", "derived from (stem_P50, stem_c), as above"},
      no_gradient{"root_b", "derived from (root_P50, root_c), as above"},
      no_gradient{"root_psi_crit", "derived from (root_P50, root_c), as above"},
      no_gradient{"a_p1", "the light-response curve the Farquhar leaf replaced"},
      no_gradient{"a_p2", "the light-response curve the Farquhar leaf replaced"},
      no_gradient{"beta1", "no equation on this path reads it"},
      no_gradient{"var_sapwood_volume_cost", "declared and carried; no equation "
                                             "on this path reads it"},
      no_gradient{"nmass_l", "declared and carried; no equation reads it"},
      no_gradient{"nmass_s", "declared and carried; no equation reads it"},
      no_gradient{"nmass_b", "declared and carried; no equation reads it"},
      no_gradient{"nmass_r", "declared and carried; no equation reads it"},
      no_gradient{"dmass_dN", "declared and carried; no equation reads it"}};

  // Every member above, in declaration order. field_ptrs() and both
  // ad_parameter_* projections read this, so a member reaches all of them or
  // none, and the sizeof assert below makes the list total.
#define PLANT_TF24_AD_PARAMETER(x) ad_parameter{&TF24_Pars::x, #x}
  static constexpr std::array ad_parameter_fields{
      PLANT_TF24_AD_PARAMETER(lma),
      PLANT_TF24_AD_PARAMETER(rho),
      PLANT_TF24_AD_PARAMETER(hmat),
      PLANT_TF24_AD_PARAMETER(omega),
      PLANT_TF24_AD_PARAMETER(eta),
      PLANT_TF24_AD_PARAMETER(theta),
      PLANT_TF24_AD_PARAMETER(a_l1),
      PLANT_TF24_AD_PARAMETER(a_l2),
      PLANT_TF24_AD_PARAMETER(a_r1),
      PLANT_TF24_AD_PARAMETER(a_cr1),
      PLANT_TF24_AD_PARAMETER(a_b1),
      PLANT_TF24_AD_PARAMETER(r_s),
      PLANT_TF24_AD_PARAMETER(r_b),
      PLANT_TF24_AD_PARAMETER(r_r),
      PLANT_TF24_AD_PARAMETER(r_cr),
      PLANT_TF24_AD_PARAMETER(r_l),
      PLANT_TF24_AD_PARAMETER(a_y),
      PLANT_TF24_AD_PARAMETER(a_bio),
      PLANT_TF24_AD_PARAMETER(k_l),
      PLANT_TF24_AD_PARAMETER(k_b),
      PLANT_TF24_AD_PARAMETER(k_s),
      PLANT_TF24_AD_PARAMETER(k_r),
      PLANT_TF24_AD_PARAMETER(k_cr),
      PLANT_TF24_AD_PARAMETER(a_p1),
      PLANT_TF24_AD_PARAMETER(a_p2),
      // Occurs only in fecundity_dt's denominator, which the fitness metrics
      // read and no census metric does: its census column is zero on any
      // trajectory. `ladder_zero_outside_the_metric_support` checks that claim.
      PLANT_TF24_AD_PARAMETER(a_f3),
      PLANT_TF24_AD_PARAMETER(a_f1),
      PLANT_TF24_AD_PARAMETER(a_f2),
      PLANT_TF24_AD_PARAMETER(S_D),
      PLANT_TF24_AD_PARAMETER(a_d0),
      PLANT_TF24_AD_PARAMETER(d_I),
      PLANT_TF24_AD_PARAMETER(a_dG1),
      PLANT_TF24_AD_PARAMETER(a_dG2),
      PLANT_TF24_AD_PARAMETER(a_st1),
      PLANT_TF24_AD_PARAMETER(a_st2),
      PLANT_TF24_AD_PARAMETER(a_st3),
      PLANT_TF24_AD_PARAMETER(storage_relaxation_offset),
      PLANT_TF24_AD_PARAMETER(a_pl0),
      PLANT_TF24_AD_PARAMETER(a_pl1),
      PLANT_TF24_AD_PARAMETER(a_pl2),
      PLANT_TF24_AD_PARAMETER(a_pl3),
      PLANT_TF24_AD_PARAMETER(a_pl4),
      PLANT_TF24_AD_PARAMETER(a_sw),
      PLANT_TF24_AD_PARAMETER(k_I),
      PLANT_TF24_AD_PARAMETER(vcmax_25),
      PLANT_TF24_AD_PARAMETER(stem_P50),
      PLANT_TF24_AD_PARAMETER(K_s),
      PLANT_TF24_AD_PARAMETER(stem_c),
      PLANT_TF24_AD_PARAMETER(stem_b),
      PLANT_TF24_AD_PARAMETER(psi_crit),
      PLANT_TF24_AD_PARAMETER(beta1),
      PLANT_TF24_AD_PARAMETER(TF24_beta2),
      PLANT_TF24_AD_PARAMETER(TF24_cost_scale),
      PLANT_TF24_AD_PARAMETER(TF24_floor_lambda_o),
      PLANT_TF24_AD_PARAMETER(jmax_25),
      PLANT_TF24_AD_PARAMETER(phi_J),
      PLANT_TF24_AD_PARAMETER(curv_fact_elec_trans),
      PLANT_TF24_AD_PARAMETER(curv_fact_colim),
      PLANT_TF24_AD_PARAMETER(R_d_25),
      PLANT_TF24_AD_PARAMETER(var_sapwood_volume_cost),
      PLANT_TF24_AD_PARAMETER(nmass_l),
      PLANT_TF24_AD_PARAMETER(nmass_s),
      PLANT_TF24_AD_PARAMETER(nmass_b),
      PLANT_TF24_AD_PARAMETER(nmass_r),
      PLANT_TF24_AD_PARAMETER(dmass_dN),
      PLANT_TF24_AD_PARAMETER(root_depth_shape_eta),
      PLANT_TF24_AD_PARAMETER(root_c),
      PLANT_TF24_AD_PARAMETER(root_P50),
      PLANT_TF24_AD_PARAMETER(root_b),
      // The root curve's dry limit, and a limit for the same reason.
      PLANT_TF24_AD_PARAMETER(root_psi_crit),
      PLANT_TF24_AD_PARAMETER(rooting_depth_max),
      // The stem path integral's two live parameters. Both reach the leaf
      // through stem_path_exponent() and the path length, so both carry a
      // column; theta_c is the third and is excluded below.
      PLANT_TF24_AD_PARAMETER(D_c),
      PLANT_TF24_AD_PARAMETER(theta_c),
      PLANT_TF24_AD_PARAMETER(L_tip),
      PLANT_TF24_AD_PARAMETER(recruitment_decay),
      PLANT_TF24_AD_PARAMETER(use_energy_balance),
      PLANT_TF24_AD_PARAMETER(d),
      PLANT_TF24_AD_PARAMETER(alpha_PAR),
      PLANT_TF24_AD_PARAMETER(alpha_NIR),
      PLANT_TF24_AD_PARAMETER(L_crown),
      PLANT_TF24_AD_PARAMETER(g_min)
  };
#undef PLANT_TF24_AD_PARAMETER

  // Every entry of the table, which is the whole parameter set: a rebind carries
  // it across a scalar change. ad_parameters() is the subset with a column.
  // Read off the two lists by plant/pars_table.h; see there for what each is.
  std::vector<S*> field_ptrs() {
    return pars_table::field_ptrs(*this, ad_parameter_fields);
  }
  static constexpr size_t field_count = ad_parameter_fields.size();
  static constexpr bool has_column(std::string_view name) {
    return pars_table::has_column(name, undifferentiable);
  }
  static constexpr size_t column_count =
      pars_table::column_count(ad_parameter_fields, undifferentiable);
  // A name matching no member excludes nothing, so the two counts stop adding
  // up: a parameter renamed or removed, or one listed twice, is a compile error
  // here rather than a list entry that has quietly stopped excluding anything.
  static_assert(column_count + undifferentiable.size() == field_count,
                "undifferentiable names a parameter ad_parameter_fields does "
                "not carry, or names one twice");
  static constexpr std::array<ad_parameter, column_count> ad_parameter_columns =
      pars_table::columns<column_count>(ad_parameter_fields, undifferentiable);
  static constexpr const auto& ad_parameter_table() { return ad_parameter_fields; }
};

// Every member of TF24_Pars is an S, so one added without extending
// ad_parameter_table() changes this size and is refused rather than dropped.
static_assert(sizeof(TF24_Pars<double>) ==
              TF24_Pars<double>::field_count * sizeof(double),
              "TF24_Pars has a member ad_parameter_fields does not list");

}

#endif
