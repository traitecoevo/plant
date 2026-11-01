// Built from  inst/include/plant/models/ff16_strategy.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_STRATEGY_H_
#define PLANT_PLANT_TF24_STRATEGY_H_

#include <plant/census.h>
#include <plant/strategy.h>
#include <odelia/ode_util.hpp>
#include <plant/models/tf24_environment.h>
#include <plant/census_gradient.h>
#include <plant/qag.h>
#include <plant/leaf_model.h>
#include <plant/canopy_shape.h>
#include <plant/stem_hydraulics.h>
#include <odelia/implicit_node.hpp>
#include <plant/with_slope.h>
#include <array>
// cstdio/cstdlib for the environment-gated curvature comparison in
// record_leaf_outputs, and nothing else here.
#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <utility>
#include <type_traits>

namespace plant {

// Biological (user-settable) parameters for the TF24 strategy. Held as a value
// member `pars` on TF24_Strategy and exposed to R as a nested RcppR6 list class
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
  // formulas -- so these two members are TF24_hyperpar's reporting copies, not
  // inputs to the model. Setting one does not change the curve.
  // stem_P50 and stem_c are TF24_hyperpar's values at the default K_s, so a
  // strategy built from traits and the bare default share one curve (#611).
  S stem_P50 = 2.8887256653360196;
  // Sapwood-specific conductivity of the TERMINAL segment. Was 1, a whole-stem
  // value under the height-linear model; back-derived to the tip by
  // TF24_K_s_from_whole_stem(1) so that resistance is UNCHANGED at
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
  // net_mass_production_dt adds `leaf.shadow_cost()` back to `leaf.profit_` before
  // converting to carbon. `Theta~` is carbon actually forgone; `lambda_o` is the
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
  // plant/stem_hydraulics.h for the closed form and
  // notes/plan-tf24-height-hydraulics.md for the derivation.
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

  // Where a parameter lives and what it is called. Nothing here says what its
  // gradient will be: the sweep computes one for every parameter except the
  // `undifferentiable` names below.
  struct ad_parameter {
    // A member pointer rather than an address in one instance, which is what
    // lets the table below be constexpr.
    S TF24_Pars::* at;
    const char* name;
  };

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
  struct no_gradient {
    const char* name;
    const char* why;
  };
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
      PLANT_TF24_AD_PARAMETER(L_crown)
  };
#undef PLANT_TF24_AD_PARAMETER

  // Every entry of the table, which is the whole parameter set: a rebind carries
  // it across a scalar change. ad_parameters() is the subset with a column.
  std::vector<S*> field_ptrs() {
    std::vector<S*> ret;
    ret.reserve(ad_parameter_fields.size());
    for (const ad_parameter& f : ad_parameter_fields) {
      ret.push_back(&(this->*f.at));
    }
    return ret;
  }

  // Read off the table, so a member added to one is added to both.
  static constexpr size_t field_count = ad_parameter_fields.size();

  // A parameter carries a column unless the list above says no gradient exists
  // for it. Derived rather than declared: adding a member gives it a column, and
  // taking one away is a deliberate entry with a sentence attached.
  static constexpr bool has_column(std::string_view name) {
    for (const no_gradient& n : undifferentiable) {
      if (name == n.name) {
        return false;
      }
    }
    return true;
  }

  // How many of them carry a gradient column. Counted here so a caller wanting
  // the width does not build a vector of pointers to measure it.
  static constexpr size_t column_count = [] {
    size_t n = 0;
    for (const ad_parameter& f : ad_parameter_fields) {
      if (has_column(f.name)) ++n;
    }
    return n;
  }();

  // A name matching no member excludes nothing, so the two counts stop adding
  // up: a parameter renamed or removed is a compile error here rather than a
  // list entry that has quietly stopped excluding anything. A name listed twice
  // fails the same way.
  static_assert(column_count + undifferentiable.size() == field_count,
                "undifferentiable names a parameter ad_parameter_fields does "
                "not carry, or names one twice");

  // The entries that carry a column, filtered once here. Filtered per call, it
  // was a string comparison against every undifferentiable name for each of the
  // sixty fields -- and the gradient asks for this list once per recording.
  static constexpr std::array<ad_parameter, column_count> ad_parameter_columns =
      [] {
        std::array<ad_parameter, column_count> ret{};
        std::size_t n = 0;
        for (const ad_parameter& f : ad_parameter_fields) {
          if (has_column(f.name)) {
            ret[n++] = f;
          }
        }
        return ret;
      }();

  static constexpr const auto& ad_parameter_table() { return ad_parameter_fields; }
};

// Every member of TF24_Pars is an S, so one added without extending
// ad_parameter_table() changes this size and is refused rather than dropped.
static_assert(sizeof(TF24_Pars<double>) ==
              TF24_Pars<double>::field_count * sizeof(double),
              "TF24_Pars has a member ad_parameter_fields does not list");

// Templated on the scalar S the state, the traits and everything derived from
// them carry; double is production. The embedded Leaf, the Control tolerances
// and the extrinsic drivers stay double.
// A cursor over the operating points one rate evaluation solves for. The points
// themselves are the SOLVER's record, beside the times, the sizes and the states,
// because they are part of what the run did; this walks the list the solver hands
// over for the evaluation about to run.
//
// At most one handle is open, and which one is the constness of what was handed
// over. Nothing here remembers a mode, so nothing can outlive the list it
// described.
// What a pass records about one leaf solve: the collar it placed, and the ARM it
// placed it on.
//
// Two values rather than a class, because that is all a recorded decision is --
// the collar MOVES and is re-derived from, the arm SELECTS and is replayed. What
// else a pass might record is rederivable from the collar alone, which is what
// phylloptim's evaluate_root_collar_psi does.
struct leaf_solved_point {
  double collar = 0.0;
  Leaf::OperatingPointKind kind = Leaf::OperatingPointKind::Unsolved;
};

class leaf_solved_points {
public:
  using points = std::vector<leaf_solved_point>;

  void store_into(points& into) {
    solved = 0;
    storing = &into;
    loading = nullptr;
    storing->clear();
  }
  void load_from(const points& from) {
    solved = 0;
    storing = nullptr;
    loading = &from;
  }
  // One evaluation is one list, so this ends where the evaluation does. A solve
  // outside a rate evaluation then stores nothing and loads nothing.
  void end() {
    solved = 0;
    storing = nullptr;
    loading = nullptr;
  }

  // The operating point the run found for the next solve of this evaluation, or one
  // holding nothing -- which is what a pass filling the record gets, and is the
  // caller's signal to solve for it.
  //
  // ⚠️ RUNNING OFF THE END IS A FAULT, not a fall back to solving. The run made one
  // entry per solve at this evaluation, so a pass that asks for more is a pass
  // whose model disagrees with the record it was handed, and searching from here
  // would tape a discretisation the run never took with every number finite.
  leaf_solved_point load() {
    if (loading == nullptr) {
      return leaf_solved_point();
    }
    if (solved >= loading->size()) {
      util::stop("leaf_solved_points: this rate evaluation asked for operating "
                 "point " + util::to_string(solved + 1) + " where the run it is "
                 "replaying solved for " + util::to_string(loading->size()));
    }
    ++placed;
    return (*loading)[solved++];
  }

  // How many points this record placed. Counted because a record that engages and
  // one that quietly does not are the same green suite: every number a placement
  // produces is the number a search produces, so nothing else can tell them apart.
  size_t placements() const { return placed; }

  void store(const leaf_solved_point& point) {
    if (storing != nullptr) {
      storing->push_back(point);
    }
  }

private:
  // The list this evaluation stores into, and the one it loads from. At most one is
  // ever open.
  points* storing = nullptr;
  const points* loading = nullptr;
  size_t solved = 0;
  size_t placed = 0;
};

// Selects TF24_Strategy's constructor from another strategy, which rebind() uses.
struct tf24_strategy_from_t {};

template <typename S = double>
class TF24_Strategy: public Strategy<TF24_Environment<S>> {
public:
  using value_type = S;

  typedef std::shared_ptr<TF24_Strategy> ptr;
  TF24_Strategy();

  // Another strategy's values at this scalar, as assign_from() writes them.
  //
  // ⚠️ THE LEAF IS COPY-CONSTRUCTED, NOT DEFAULTED AND THEN ASSIGNED. A default
  // phylloptim::Leaf computes both vulnerability tables from incomplete gamma
  // functions, about 60 times the cost of copying them, and assign_from() would
  // overwrite them at once. The sweep rebinds the patch on the order of once
  // per recorded step, so each build it skips is one per step.
  template <class S1>
  TF24_Strategy(tf24_strategy_from_t, const TF24_Strategy<S1>& src)
      : leaf(src.leaf) {
    assign_from(src);
  }

  // Scientific version. Bump ONLY when equations or default parameters change
  // the simulation output for identical inputs. Do NOT bump for refactors,
  // performance, interface, or serialisation changes. Bumping invalidates
  // logpile's cache for this model (see plant::model_version() / model_id()).
  // Starts at 2: a published result exists using pre-versioning "v1" science.
  // v3 (#517/#554): NSC storage pool + reserve-gated growth and reserves-based
  // mortality change the simulation output for identical inputs; TF24f's
  // compound version auto-tracks this to 3.1.
  // v4: the hydraulic shut-down exits no longer leave the previous step's
  // transport state in place, so a shut-down plant stops drawing water from the
  // soil (it used to keep extracting its last wet-step uptake, because `Leaf`
  // is reused across steps and soil_consumption_ feeds the patch water
  // balance). Water-limited runs therefore change: on the scenario gateway,
  // offspring production moves by up to 5e-3 relative on 5 of 8 scenarios,
  // while every success/failure classification is unchanged. TF24f's compound
  // version auto-tracks this to 4.1.
  // v5: the leaf gas-exchange and hydraulics model is now the standalone
  // `phylloptim` package rather than a copy in this repo, and the swap carries four science
  // changes. Measured on the one-species SCM scenario of test-strategy-tf24.R
  // (max_patch_lifetime = 5), offspring production moves 81.9083 -> 83.9026,
  // i.e. **+2.4%**. Attributed by re-running both arms with the atm_kpa driver
  // forced to 101.3, which removes the pressure change and leaves the rest:
  //
  //     arm                  atm_kpa 100.5     atm_kpa 101.3
  //     this repo's leaf        81.9083           81.8201
  //     the leaf package        83.9026           81.8985
  //
  // So **the pressure fix is ~25x the rest of the swap put together** (+2.4%
  // against +0.10%), which was not the expectation going in. The leaf package's
  // ppm-to-Pa conversion is derived from atm_kpa (phylloptim #15 item 10c) instead
  // of hard-coded at 0.1013 = 1e-6 * 101300 Pa; TF24_Environment's atm_kpa driver
  // defaults to **100.5**, so Gamma*, Kc, Ko, Km and the ci root-find bounds all
  // move. Before the fix the conductance side of the model responded to atm_kpa
  // while the photosynthesis side silently assumed sea level -- visible in the
  // table as this repo's leaf moving only -0.11% across the same 0.8 kPa that
  // moves the package -2.4%.
  //
  // The remaining +0.10% is the two further stale-state exits (phylloptim #26,
  // ported from #585) plus the supply-path extraction (phylloptim #2). TF24f's
  // compound version auto-tracks this to 5.1.
  // v6: the leaf package moved to ONE representation for water potential --
  // positive magnitudes throughout (phylloptim #25). Two consequences, and the first
  // is why this is a version bump rather than a refactor:
  //   * the **`opt_root_psi` aux changes sign**. It is now the positive magnitude,
  //     which is what TF24f's `opt_root_psi_state` has always held. Before, the aux
  //     reported the signed potential while tf24f_strategy.cpp negated it back for
  //     the state -- an inconsistency in plant's own reported outputs, and the two
  //     compensating negations are deleted here. Any stored output or cached
  //     analysis reading that aux would silently change meaning, which is exactly
  //     what model_version() exists to catch.
  //   * outputs move slightly: one-species SCM offspring production 83.9026 ->
  //     83.8761, i.e. **-3.2e-4 relative**. This is NOT an equation change. The
  //     rewrite is exactly sign-symmetric in IEEE; what is not is boost's TOMS748,
  //     whose iterates depend on the bracket's orientation, and #25 reverses it.
  //     Measured there: 12 of 288 golden operating points differ by 1-3 ULP, the
  //     rest exactly. The SCM's adaptive stepper and node schedule amplify that to
  //     3e-4 -- within the ~GSS_tol_abs (1e-3) ceiling the leaf package documents,
  //     and small enough that every pinned test value and the exact stochastic
  //     TF24 counts (101/23) pass unchanged.
  // TF24f's compound version auto-tracks this to 6.1.
  // v7: the collar bracket is finally clamped to root_psi_crit, the potential at
  // which root conductivity is down to 5% (phylloptim #24, plant #584). The clamp was
  // written as a std::max against a *signed* root_psi_crit, so it could never bind
  // and the solver optimised over a collar the root system cannot supply. **The
  // window is 1.2 MPa wide at TF24's defaults** -- psi_crit = 7.085493 against
  // root_psi_crit = 5.870283 -- so this is a dry-corner correction, not a rounding
  // one. Two regimes inside it: the interval is tightened (the plant still
  // transpires, at a wetter collar), or root_psi_crit lands below the zero-uptake
  // collar and the plant shuts down because no operating point both moves water and
  // stays inside the root limit.
  //
  // No tested scenario moves: the standard SCM run stays at 83.8761 offspring, bit
  // for bit, because a mesic patch never drives the collar past 5.87 MPa. It is
  // still a version bump -- a user running a dry scenario gets different (correct)
  // numbers for identical inputs, and logpile's cache has to know.
  // TF24f's compound version auto-tracks this to 7.1.
  // v8: `TF24_Environment`'s `atm_kpa` driver default goes **100.5 -> 101.3**, and
  // this is the entry to read if you only read one. It largely CANCELS v5.
  //
  // v5 recorded +2.4% from deriving the leaf's ppm -> Pa conversion from `atm_kpa`
  // instead of hard-coding 0.1013. That constant *was* 101.3 kPa in disguise
  // (1e-6 * 101300 Pa), so the shift was not the fix doing damage -- it was this
  // driver disagreeing with the rest of the model. 100.5 arrived in `34d46ac2`
  // ("Simplify scm & environment interface", #446), an interface refactor that does
  // not mention atmospheric pressure, with no rationale recorded anywhere, while
  // every leaf-level test used 101.3. An artefact, not a site elevation.
  //
  // Pinning it to the value the model already assumed collapses the whole branch's
  // movement. Net effect of ALL of it (the swap, the #15 catch-up, the #26 ported
  // fixes, #25 and #24) against `develop`:
  //
  //     one-species SCM offspring   81.9083 -> 81.7426     -0.20%
  //     stochastic TF24 counts      103 / 28 -> 103 / 28   unchanged, exactly
  //
  // So **+2.43% became -0.20%**, and every pinned baseline reverts to develop's own
  // values: the two SCM offspring figures pass at their original 82.09077702 /
  // 67.54060383, and the seeded stochastic integers match bit for bit. That exact
  // match on discrete counts is a sharper statement than any tolerance-based check
  // that the swap preserves TF24's science.
  //
  // The fix itself is NOT undone -- an off-sea-level run still gets a self-consistent
  // Gamma*/Kc/Ko/Km and conductance side, which is the whole point of item 10c. Set
  // `atm_kpa` per site if you mean altitude; it just no longer defaults to an
  // altitude nobody chose.
  // v9: the storage pool's rate is a charge and a drain limited separately, so
  // the pool stays within capacity by the shape of its own flow. Storage
  // previously had no upper bound at all -- the read was clipped at capacity
  // while the state ran to 1.035 of it, with half a full-lifetime stand sitting
  // at or above the clip -- and that surplus now stays in production instead.
  // v9 also: the light field and the census take their trapezium widths from the
  // coordinate the density is carried in. TF24 defaults to the birth-date
  // coordinate, so this moves its output; the water reduction already integrated
  // over birth dates and is unchanged.
  //
  // v10 (#615): stem resistance becomes a path integral over two within-plant
  // anatomical profiles instead of being linear in height, so the height
  // exponent is derived rather than assumed. Defaults D_c = 0.2, L_tip = 0.02,
  // theta_c = 0, giving R_L ~ H^0.6. K_s is reparameterised from the old
  // whole-stem 1 to the terminal-segment 0.33577377016801868 (a factor 2.9782)
  // so that resistance is UNCHANGED at TF24_H_ANCHOR = 1 m. theta_c stays at 0,
  // so `theta` keeps its whole-plant meaning and no parameter file needs
  // migrating.
  //
  // Only theta/K_s is identifiable, so there is one free scalar and the two
  // models agree at exactly ONE height. It is a ROTATION about the anchor:
  //
  //     H (m)        0.394  1.00   5.00   8.00   16.60  30.0   60.0
  //     R_new/R_old  1.38   1.00   0.55   0.46   0.35   0.28   0.21
  //
  // Only sub-metre plants pay more than they did; everything taller pays
  // progressively less. Setting D_c, theta_c and L_tip to zero and K_s to 1
  // recovers the previous model exactly, end-to-end through the SCM -- see the
  // "height-linear parameters" test in tests/testthat/test-strategy-tf24.R,
  // which reproduces the pinned values from before this change, unmodified.
  //
  // THE SIGN OF THE EFFECT DEPENDS ON DENSITY, which is the most important
  // thing to know about this change. Individually, plants are better off
  // wherever they are taller than the anchor (assimilation ratio 0.9968 at
  // 0.5 m, 1.0006 at 1 m, 1.0543 at 5 m, 1.1410 at 10 m, single plant, wet soil,
  // no competition). But lower resistance also means faster transpiration, so in a
  // dense stand everyone draws the shared soil column down faster and the patch
  // does WORSE. One-species SCM, hmat = 5, max_patch_lifetime = 5:
  //
  //     birth_rate    0.5      2       20
  //     ratio        1.196   1.047   0.803
  //
  // The pinned scenarios below all run at birth_rate = 20, i.e. at the least
  // favourable end of that range; they are not representative of the change's
  // sign in general.
  //
  //     one-species SCM offspring      30.2980 ->  24.3214   -19.73%
  //     two-species, fast              23.2557 ->  18.5430   -20.26%
  //     two-species, slow            4.0284e-6 -> 1.8965e-6  -52.92%
  //     birth-date coordinate, fast   233.3660 -> 219.2668    -6.04%
  //     birth-date coordinate, slow    43.7240 ->  34.5877   -20.90%
  //     seeded stochastic counts         79 / 3 ->  77 / 3
  //
  // The hydraulic gateway runs longer patches at the default hmat and moves the
  // other way, up by 3.1x to 2208x on every scenario, with S01 and S02 crossing
  // R0 = 1 so persistence goes 1/8 -> 3/8. 8/8 still run, 0 crash.
  //
  // theta_c is declared but REFUSED (prepare_strategy throws on any non-zero
  // value). theta is read by the carbon budget as well as the hydraulic term,
  // so a hydraulics-only profile would be an incoherent model rather than a
  // staging step; it lands everywhere at once or not at all.
  //
  // stem_P50 is deliberately UNCHANGED (2.8887 MPa): make_TF24_hyperpar derives
  // the vulnerability curve from K_s, so B_Hv1 was re-anchored 0.4607063 ->
  // 0.36591565341924093 to stop the reparameterisation from also moving it.
  // Left alone it would have gone to 3.5933 MPa.
  // v11: three changes that move output at identical inputs.
  //   * R_d_25 = 1.44, leaf dark respiration at 25 C, is a NEW parameter; it did
  //     not exist on develop and it reaches phylloptim as par_R_d_25.
  //   * The root vulnerability curve is reparameterised to (P50, c), matching
  //     the stem side and phylloptim, so root_P50 = 3.4 is the free parameter
  //     and root_b is derived from it. root_b 3.898245 -> 3.8982451221145307
  //     and root_psi_crit 5.8702825428827037 -> 5.8702827267723245, both
  //     +3.1e-10 relative.
  //   * vulnerability_curve_ncontrol 100 -> 400, tracking phylloptim's
  //     Leaf::ncontrol_default, so the pre-tabulated curve the root-finds run
  //     on is four times finer.
  //
  // FF16 and K93 are NOT bumped. Their snapshots move on the same
  // vulnerability_curve_ncontrol and on the new gradient_curvature_floor, and
  // neither model reads either: the curve count is read at one site, building
  // the Leaf, and the floor only by the sweep. Both suites are green on their
  // pinned outputs.
  // v12: as FF16 v3 -- odelia 0.8.0's step-size rule and the default ode_tol
  // 1e-4 -> 1e-5 move every adaptive run within the ODE tolerance. TF24f
  // follows to v12.1.
  // v13: the light chain. The driver is SWdown (W m^-2, default 900) where it
  //   was PPFD (1800), so the default PPFD above the canopy is 2056.5; the
  //   leaf absorbs its intercepted light itself (alpha_PAR 0.85, alpha_NIR
  //   0.20), so on the energy-balance path net radiation falls to about half;
  //   and `a` is `phi_J` on an absorbed basis, which reassociates J. TF24f and
  //   TF24t follow to v13.1 and v13.4.
  // v14: phylloptim 0.16.0 -> 0.26.2, the leaf where it stops moving water or
  //   gas. A shut leaf is the ordinary evaluator at its wet bound, so in light it
  //   nets zero instead of -R_d and profit no longer jumps at shutdown or at the
  //   light-compensation point; psi_crit is no longer a bound, so the stem
  //   follows a drying soil past it and the stands' dry pins at psi_crit become
  //   interior points; shade death is gone as a kind; the vulnerability integral
  //   is a C2 table; and the ci root-find is solved to 1e-13. Output moves near
  //   the root limit, below light compensation and in heat. TF24f, TF24t and
  //   TF24w follow; TF24t's thermal layer also moves (thermal_revision 5).
  static constexpr int scientific_version = 14;

  // calculate the amount of water transpired relativised by leaf area index.

  S evapotranspiration_dt(const S& area_leaf_, int soil_layer);


  // Overrides ----------------------------------------------

  // update this when the length of state_names changes
  static size_t state_size () { return 9; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
      "area_heartwood",
      "mass_heartwood",
      "storage",
      "log_area_leaf_departure",
      "log_area_sapwood_departure",
      "leaf_marginal_return_tracked"
      });
  }

  // TF24 carries reserves that feed back into growth, so its growth rate is not
  // a function of size alone and the height coordinate's compression term is a
  // different derivative rather than a worse approximation of the right one
  // (#590). The two coordinates' offspring production differ by several times
  // here, where FF16's and K93's agree to 5e-5 and 5.6e-7, so this is a
  // correctness default and not a performance one. TF24f, TF24t and TF24w
  // inherit it.
  static constexpr bool density_in_birth_date_default = true;

  // Each species' net reproduction ratio and offspring production are stand
  // metrics (SCM::stand_metrics), refereed in test-tf24-fitness-gradient.R.
  static constexpr bool differentiable_fitness = true;

  // The metrics a census of this model sums, in the order it reports them. Beside
  // the state and aux names because all three say what this model calls its own
  // quantities -- and every slot is read by its cached index, which is how every
  // other reader of these slots reaches them.
  //
  // A table rather than a list built per call, for the reason
  // ad_parameter_fields is one: a census is read from several entry points and
  // each was allocating its own copy of this.
  //
  // ⚠️ THE LEAF AREA IS THE ONE THE PLANT HAS, read from the departure state, not
  // the one its height implies, and the sapwood is likewise at the Huber value it
  // carries (#516). At rest both departures are exactly zero and exp(0) is
  // exactly one, so these are the height-only expressions bit for bit.
  static const auto& census_metrics() {
    using strategy = TF24_Strategy<S>;
    static constexpr std::array<census_metric<strategy>, 3> metrics{{
      {"leaf_area",
       [](const strategy& p, const Internals<S>& vars) -> S {
         return p.area_leaf_carried(vars);
       }},
      {"mass_above_ground",
       [](const strategy& p, const Internals<S>& vars) -> S {
         const S height = vars.state(HEIGHT_INDEX);
         const S area_leaf = p.area_leaf_carried(vars);
         return p.mass_above_ground(
             p.mass_leaf(area_leaf),
             p.mass_bark(p.area_bark(area_leaf), height),
             p.mass_sapwood(
                 p.area_sapwood(
                     area_leaf,
                     vars.state(p.state_idx_log_area_sapwood_departure)),
                 height),
             vars.state(p.state_idx_mass_heartwood));
       }},
      {"area_stem",
       [](const strategy& p, const Internals<S>& vars) -> S {
         const S area_leaf = p.area_leaf_carried(vars);
         return p.area_stem(
             p.area_bark(area_leaf),
             p.area_sapwood(area_leaf,
                            vars.state(p.state_idx_log_area_sapwood_departure)),
             vars.state(p.state_idx_area_heartwood));
       }},
    }};
    return metrics;
  }

  // The leaf area a plant carries: A*(h) exp(phi). update_dependent_aux keeps the
  // same number in the competition_effect aux; this reads the states directly,
  // for a caller holding states whose aux may not have been refreshed.
  S area_leaf_carried(const Internals<S>& vars) const {
    using std::exp;
    return area_leaf(vars.state(HEIGHT_INDEX)) *
      exp(vars.state(state_idx_log_area_leaf_departure));
  }

  // The storage pool. Its rate holds the exact flow inside [0, S_max], so a
  // negative value is a step that overshot rather than a state the model has,
  // and the solver rejects it and retries smaller.
  static std::vector<std::string> non_negative_states() { return {"storage"}; }

  // Virtual so that a subclass can append its own readings (TF24w's per-column
  // fluxes, #681): refresh_indices builds aux_index from this, and aux_size()
  // counts it, so one override reaches both.
  virtual std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "height_inverse",
      "net_mass_production_dt",
      "root_mass",
      "opt_psi_stem",
      "opt_root_psi",
      "transpiration",
      "E_up_",
      "profit",
      // The part of the cost the objective deducted that is a PRICE rather than
      // a realised carbon loss (umol C m^-2 s^-1): `TF24_floor_lambda_o * E` on
      // the seated curve, and exactly zero at the default price. Reported beside
      // `profit` because without it the two cannot be separated after the fact
      // from plant's own output, and the distinction is the whole point of the
      // curve: the carbon the plant actually kept is `profit + shadow_cost`,
      // which is what net_mass_production_dt grows on.
      "shadow_cost",
      "stom_cond_CO2",
      // Net CO2 assimilation at the optimal operating point, per unit leaf
      // area (umol CO2 m^-2 s^-1). Net, not gross: Leaf::assim_colimited()
      // subtracts dark respiration R_d_, so gross = assimilation + R_d_ with
      // R_d_ = 0.015 * vcmax_ at the acclimated vcmax_.
      "assimilation",
      // Leaf temperature at the optimal operating point (deg C) -- an OUTPUT,
      // not the `leaf_temp` driver, and the distinction is the point of
      // reporting it. With pars.use_energy_balance off this equals the driver,
      // so the column is flat at TF24's defaults; with it on the leaf solves
      // its own temperature from its transpiration per operating point, and
      // that value was previously computed, used to re-derive the whole
      // Farquhar block, and then discarded -- so the one quantity the
      // Penman-Monteith path exists to produce was the one a canopy-level
      // analysis could not read (#625). Reported here rather than left to be
      // inferred from an assimilation that cannot be explained without it.
      //
      // Under the deep-crown shading model this is the leaf-area-weighted crown
      // mean, integrated alongside the other leaf outputs; it is NOT the
      // temperature of any single leaf, and a canopy with a hot top and a cool
      // base reports the mean of the two. The depth profile itself is not an
      // aux (a fixed-width slot cannot carry a per-quadrature-node vector).
      "Tleaf",
      // The marginal leaf's carbon balance, which the replacement gate reads.
      // Dimensionless and zero at break-even, so its sign says whether the last
      // leaf is paying for itself.
      "leaf_marginal_return",
      // The marginal return of SAPWOOD area, in growth rather than production
      // (#516). Dimensionless, and zero at the Huber value that maximises the
      // plant's growth rate, so its sign says which way the stem should move.
      // Zero unless sapwood acclimation is on (a_sw > 0), because computing it
      // costs two further leaf solves.
      "sapwood_marginal_return"
    });
    // add the associated computation to compute_rates and compute there
    if (this->collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  using ad_parameter = typename TF24_Pars<S>::ad_parameter;

  // A constant birth rate, at the working scalar so that it is a gradient
  // column; a variable one is a driver read as a double and has no column.
  S birth_rate_ = 1.0;
  bool has_birth_rate_column() const { return !this->is_variable_birth_rate; }
  S birth_rate_at(double t) const {
    return has_birth_rate_column()
             ? birth_rate_
             : S(this->extrinsic_drivers.evaluate("birth_rate", t));
  }

  // How many gradient columns one of these carries, for a caller sizing a batch.
  size_t ad_column_count() const {
    return TF24_Pars<S>::column_count + (has_birth_rate_column() ? 1 : 0);
  }

  // Addresses of the parameters a gradient can be taken with respect to, in the
  // order ad_parameter_names() gives: the table entries whose role has a column,
  // then the birth rate where it is one. Both allocate, so take them once per
  // gradient evaluation and hold them for the run rather than per block; the
  // strategy is shared and the fields do not move. Index against .size().
  std::vector<S*> ad_parameters() {
    std::vector<S*> ret;
    ret.reserve(ad_column_count());
    for (const ad_parameter& p : TF24_Pars<S>::ad_parameter_columns) {
      ret.push_back(&(pars.*p.at));
    }
    if (has_birth_rate_column()) {
      ret.push_back(&birth_rate_);
    }
    return ret;
  }

  // The name each gradient column carries, one per ad_parameters() entry and in
  // that order. Which parameters are left out, and why, is `undifferentiable`.
  std::vector<std::string> ad_parameter_names() {
    std::vector<std::string> ret;
    ret.reserve(ad_column_count());
    for (const ad_parameter& p : TF24_Pars<S>::ad_parameter_columns) {
      ret.push_back(p.name);
    }
    if (has_birth_rate_column()) {
      ret.push_back("birth_rate");
    }
    return ret;
  }

  // The parameters no gradient exists for, each with the sentence that says
  // why, so a caller asking for one is refused in those words rather than told
  // the name is unknown.
  static std::vector<std::pair<std::string, std::string>> undifferentiable_reasons() {
    std::vector<std::pair<std::string, std::string>> ret;
    for (const auto& n : TF24_Pars<S>::undifferentiable) {
      ret.emplace_back(n.name, n.why);
    }
    return ret;
  }

  // Translate generic methods to TF24 strategy leaf area methods

  S competition_effect(const S& height) const {
    return area_leaf(height);
  }

  void refresh_indices();


  // TF24 Methods  ----------------------------------------------

  // [eqn 2] area_leaf (inverse of [eqn 3])
  S area_leaf(const S& height) const;

  // [eqn 1] mass_leaf (inverse of [eqn 2])
  S mass_leaf(const S& area_leaf) const;

  // Sapwood area per unit leaf area -- the Huber value, `theta` displaced by
  // the sapwood departure state. ALWAYS form the ratio this way rather than as
  // area_sapwood / area_leaf: the two are algebraically equal but not bitwise,
  // because a multiply followed by a divide does not recover its operand, and
  // that difference is enough to lose the exactness of the resting model.
  S sapwood_per_leaf_area(const S& sapwood_departure) const {
    using std::exp;
    return pars.theta * exp(sapwood_departure);
  }

  // [eqn 4] area and mass of sapwood. The one-argument form is the resting
  // allometry (used where a plant is on its trajectory by definition, e.g. the
  // germination inversion); the two-argument form reads the departure state.
  S area_sapwood(const S& area_leaf) const;
  S area_sapwood(const S& area_leaf, const S& sapwood_departure) const;
  S mass_sapwood(const S& area_sapwood, const S& height) const;

  // [eqn 5] area and mass of bark
  S area_bark(const S& area_leaf) const;
  S mass_bark (const S& area_bark, const S& height) const;

  S area_stem(const S& area_bark, const S& area_sapwood,
                            const S& area_heartwood) const;
  S diameter_stem(const S& area_stem) const;

  // [eqn 7] Mass of (fine) roots
  S mass_root(const S& area_leaf) const;

  // [eqn 7b] Mass of coarse (structural) roots, as a fraction of sapwood mass
  // (see TF24_Pars::a_cr1). Zero unless a_cr1 is set.
  S mass_coarse_root(const S& mass_sapwood) const;

  // [eqn 8] Total Mass
  S mass_live(const S& mass_leaf, const S& mass_bark,
              const S& mass_sapwood, const S& mass_root,
              const S& mass_coarse_root) const;

  S mass_total(const S& mass_leaf, const S& mass_bark, const S& mass_sapwood,
               const S& mass_heartwood, const S& mass_root,
               const S& mass_coarse_root) const;

  // Above-ground mass = leaf + all stem components (bark + sapwood +
  // heartwood); excludes roots.
  S mass_above_ground(const S& mass_leaf, const S& mass_bark,
                      const S& mass_sapwood, const S& mass_heartwood) const;

  void compute_rates(const TF24_Environment<S>& environment,
                Internals<S>& vars);
  
  void compute_roots(const TF24_Environment<S>& environment,
                Internals<S>& vars);

  void update_dependent_aux(const int index, Internals<S>& vars);

  // * Mass production
  // [eqn 13] Total maintenance respiration
  S respiration(const S& mass_leaf, const S& mass_sapwood,
                const S& mass_bark, const S& mass_root,
                const S& mass_coarse_root) const;

  S respiration_leaf(const S& mass) const;
  S respiration_bark(const S& mass) const;
  S respiration_sapwood(const S& mass) const;
  S respiration_root(const S& mass) const;
  S respiration_coarse_root(const S& mass) const;

  // [eqn 14] Total turnover
  // Turnover the plant actually PAYS for. `replacement_leaf` and
  // `replacement_sapwood` are the fractions of the leaf and sapwood flux it
  // chooses to rebuild; the rest is not charged, and the matching tissue is lost
  // through the departure states instead.
  //
  // ⚠️ ONLY A POOL WHOSE SIZE IS A STATE MAY BE GATED. Bark and fine-root
  // turnover stay fully charged because neither has a state to lose, so
  // withholding their cost would save carbon that no tissue paid for -- a free
  // lunch, and a silent one. Coarse roots follow the sapwood's mass, so they do
  // shrink with it, but they stay charged in full: over-charging is the safe
  // side, and the pool is empty at the default a_cr1 = 0.
  S turnover(const S& replacement_leaf, const S& replacement_sapwood,
             const S& mass_leaf, const S& mass_bark,
             const S& mass_sapwood, const S& mass_root,
             const S& mass_coarse_root) const;
  S turnover_leaf(const S& mass) const;
  S turnover_bark(const S& mass) const;
  S turnover_sapwood(const S& mass) const;
  S turnover_root(const S& mass) const;
  S turnover_coarse_root(const S& mass) const;

  // The fraction of leaf and sapwood turnover a plant rebuilds, as a smooth
  // function of the tracked marginal leaf return: ~1 while the marginal leaf
  // pays for itself, falling to 1 - a_pl0 as it stops. a_pl0 = 0 gives exactly 1
  // and so exactly the fixed allometry.
  //
  // Two forms of one logistic, chosen by the sign of the exponent, because a
  // healthy plant's return sits tens of gate widths above the centre: there
  // e^x overflows to inf, which the value survives and its derivative,
  // inf/inf, does not -- and a NaN partial multiplied by a_pl0 = 0 is still NaN.
  S replacement_fraction(const S& marginal_return) const {
    using std::exp;
    const S x = (marginal_return - pars.a_pl1) / plasticity_gate_width;
    if (x > 0.0) {
      const S e = exp(-x);
      return 1.0 - pars.a_pl0 * e / (1.0 + e);
    }
    return 1.0 - pars.a_pl0 / (1.0 + exp(x));
  }

  // How much of the withheld replacement actually happens, tapering to nothing
  // at the canopy floor. Same functional form as the rebuild share, mirrored:
  // rebuilding vanishes at the ceiling (phi = 0) and shedding at the floor
  // (phi = ln a_pl3), so phi is confined to [ln a_pl3, 0] by the SHAPE of the
  // flow rather than by any clamp -- and both bounds are therefore differentiable
  // and cost the solver nothing.
  S shedding_gate(const S& leaf_departure) const {
    using std::exp;
    using std::log;
    const S floor_departure = log(pars.a_pl3);
    return 1.0 - exp((floor_departure - leaf_departure) / pars.a_pl2);
  }

  // Sapwood replaces a scaled share of what its conversion to heartwood
  // removes: the same carbon gate as leaves, times exp(-psi / a_pl2), so a stem
  // carrying more conducting area than its canopy needs declines to renew the
  // excess even when carbon is ample.
  //
  // This is what gives the Huber value a restoring force. Without it psi has
  // only terms that push it around and none that return it, and measured, it
  // drifted to leave plants hydraulically UNDER-built (#516).
  //
  // ⚠️ psi >= 0 IS INVARIANT, and that is the reason this form needs no bound.
  // At psi = 0 the factor is 1, so the two replacement fractions coincide and
  //     dpsi/dt = (1 - replacement) * (k_l - k_s) >= 0,
  // since k_l > k_s. The flow at the boundary points inward, so the "too little
  // sapwood" branch is unreachable, the factor stays within (0, 1], and no
  // over-replacement (a fraction above 1) can arise to cost unbounded carbon.
  // If k_s were ever raised above k_l that argument reverses and this needs
  // revisiting.
  S sapwood_replacement_fraction(const S& replacement_leaf,
                                 const S& sapwood_departure) const {
    using std::exp;
    return replacement_leaf * exp(-sapwood_departure / pars.a_pl2);
  }

  // Carbon cost of a unit of leaf area bought at UNCHANGED height: the leaf,
  // the fine root, and the sapwood and bark that supply it, at the ratio the
  // plant prefers rather than the one it currently has. Buying the whole
  // package is what stops rebuilding a canopy from diluting the Huber value.
  S mass_live_per_area_leaf_at_height(const S& height) const {
    return pars.lma + pars.a_r1 +
      pars.theta * height * eta_c * pars.rho * (1.0 + pars.a_b1);
  }

  // [eqn 15] Net production
  S net_mass_production_dt_A(const S& assimilation, const S& respiration,
                             const S& turnover) const;

  // Takes the three #516 states beside the sizes: the sapwood departure moves
  // the conductance handed to the leaf and the sapwood carried, and the leaf
  // departure and the tracked marginal return set what turnover is charged. A
  // plant on its allometry with no history passes zeros.
  virtual S net_mass_production_dt(const TF24_Environment<S>& environment,
                                const S& height, const S& area_leaf_,
                                const S& height_inverse,
                                const S& sapwood_departure,
                                const S& leaf_departure,
                                const S& tracked_marginal);

  // The soil water potentials the leaf's supply solve reads, one per layer.
  // TF24 reads the patch's shared column; TF24w (#681) overrides this to hand
  // over the cohort's own root-zone column, computed from its own water states
  // before the base compute_rates runs. Virtual for the reason solve_leaf is:
  // net_mass_production_dt is reused unchanged by the subclass.
  virtual const std::vector<S>& soil_water_potential_for_leaf(
      const TF24_Environment<S>& environment) {
    return environment.get_soil_water_potential_state();
  }

  // Resolve the leaf operating point on the already-set-up `leaf` (i.e. after
  // leaf.set_physiology(...)). Base TF24 optimises the root-collar psi via
  // golden-section search; the TF24f variant overrides this to make the optimum
  // chase a tracked ODE state (#525). Called per crown light point from
  // net_mass_production_dt, so it must be virtual to dispatch to the override
  // when net_mass_production_dt is reused unchanged by the subclass.
  virtual void solve_leaf();

  // Read how the leaf's outputs respond to what it was given, and record each
  // output that re-enters the active chain carrying that response. The leaf is
  // already driven and solved; nothing here re-supplies it.
  void record_leaf_outputs(const S& PPFD_int, const std::vector<S>& psi_soil,
                           const S& conductance_max);
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // height state and the cached aux slots itself, so the generic Individual
  // does not need to know TF24's state/aux layout.
  S net_mass_production_dt(const TF24_Environment<S>& environment,
                                const Internals<S>& vars) {
    return net_mass_production_dt(
        environment, vars.state(HEIGHT_INDEX),
        vars.aux(aux_idx_competition_effect),
        vars.aux(aux_idx_height_inverse),
        vars.state(state_idx_log_area_sapwood_departure),
        vars.state(state_idx_log_area_leaf_departure),
        vars.state(state_idx_leaf_marginal_tracked));
  }

  // [eqn 16] Fraction of whole plan growth that is leaf
  virtual S fraction_allocation_reproduction(const S& height) const;

  // The height above which this strategy allocates at least 1 - epsilon of its
  // maximum, a_f1, to reproduction: where a continuous patch's absorbing top
  // layer sits (#519). Growth's share of surplus there is 1 - (1 - epsilon) a_f1,
  // within epsilon a_f1 of its floor -- and epsilon itself only where a_f1 is one.
  // The inverse of fraction_allocation_reproduction()'s logistic, and NaN where
  // the logistic does not rise with height (a_f2 <= 0) or hmat is not positive.
  double asymptotic_height(double epsilon) const {
    using odelia::util::to_passive;
    const double hmat = to_passive(pars.hmat);
    const double a_f2 = to_passive(pars.a_f2);
    if (!(hmat > 0.0) || !(a_f2 > 0.0)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return hmat * (1.0 + std::log((1.0 - epsilon) / epsilon) / a_f2);
  }
  S fraction_allocation_growth(const S& height) const;
  // [eqn 17] Rate of offspring production
  S fecundity_dt(const S& net_mass_production_dt,
                 const S& fraction_allocation_reproduction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  // Leaf area bought per unit live mass. Takes the sapwood departure because
  // the stem that comes with a new leaf is priced at the Huber value the plant
  // ACTUALLY carries, not at pars.theta. Without that a plant holding extra
  // conducting area gets its hydraulic benefit free, which is most of why net
  // production appeared to peak at 3.3x the pipe-model ratio.
  S darea_leaf_dmass_live(const S& area_leaf,
                          const S& sapwood_departure) const;

  // change in height per change in leaf area
  S dheight_darea_leaf(const S& area_leaf) const;
  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  S dmass_leaf_darea_leaf(const S& area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  S dmass_sapwood_darea_leaf(const S& area_leaf,
                             const S& sapwood_departure) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  // Bark is pinned to LEAF area, not to sapwood, so it keeps pars.theta and
  // does not move with the departure. Sharing the sapwood expression would give
  // it the Huber value's excursion, which it has no business having.
  S dmass_bark_darea_leaf(const S& area_leaf) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  S dmass_root_darea_leaf(const S& area_leaf) const;
  // Mass of coarse root needed for new unit area leaf, d m_cr / d a_l.
  // == a_cr1 * dmass_sapwood_darea_leaf, so it carries the same height scaling,
  // and the same Huber value: coarse roots continue the sapwood the plant has.
  S dmass_coarse_root_darea_leaf(const S& area_leaf,
                                 const S& sapwood_departure) const;
  // Growth rate of basal diameter_stem per unit stem area
  S ddiameter_stem_darea_stem(const S& area_stem) const;
  // Growth rate of components per unit time:
  S area_leaf_dt(const S& area_leaf_dt) const;
  S area_sapwood_dt(const S& area_leaf_dt) const;
  S area_heartwood_dt(const S& area_leaf, const S& sapwood_departure) const;
  S area_bark_dt(const S& area_leaf_dt) const;
  S area_stem_dt(const S& area_leaf, const S& area_leaf_dt) const;
  S diameter_stem_dt(const S& area_stem, const S& area_stem_dt) const;
  S mass_root_dt(const S& area_leaf,
                 const S& area_leaf_dt) const;
  S mass_coarse_root_dt(const S& area_leaf,
                        const S& area_leaf_dt) const;
  S mass_live_dt(const S& fraction_allocation_reproduction,
                 const S& net_mass_production_dt) const;
  S mass_total_dt(const S& fraction_allocation_reproduction,
                  const S& net_mass_production_dt,
                  const S& mass_heartwood_dt) const;
  S mass_above_ground_dt(const S& area_leaf,
                         const S& fraction_allocation_reproduction,
                         const S& net_mass_production_dt,
                         const S& mass_heartwood_dt,
                         const S& area_leaf_dt) const;

  S mass_heartwood_dt(const S& mass_sapwood) const;

  S mass_live_given_height(const S& height) const;
  S height_given_mass_leaf(const S& mass_leaf_) const;


  S mortality_dt(const S& relative_reserves, const S& cumulative_mortality) const;
  S mortality_growth_independent_dt()const ;
  // Storage-dependent growth mortality (#517): rises smoothly as relative
  // reserves r = S/S_max deplete, bounded in [a_dG1*e^-a_dG2, a_dG1].
  S mortality_storage_dependent_dt(const S& relative_reserves) const;
  // NSC storage capacity S_max = a_st1 * mass_sapwood [kg NSC], at the
  // sapwood the plant carries.
  S storage_capacity(const S& area_leaf, const S& height,
                     const S& sapwood_departure) const;
  // Initialise the storage state for a newly germinated individual (#517).
  void set_initial_states(const TF24_Environment<S>& environment, Internals<S>& vars);
  // [eqn 20] Survival of seedlings during establishment, from the carbon a
  // seedling produces at birth size. This form works that carbon out.
  S establishment_probability(const TF24_Environment<S>& environment);
  // The same, for a newborn whose rates have just been computed. A newborn is
  // already at birth size, so compute_rates has left that carbon in aux and the
  // leaf need not be solved there twice.
  S establishment_probability(const TF24_Environment<S>& environment,
                              const Internals<S>& vars) {
    return establishment_probability(environment,
                                     vars.aux(aux_idx_net_mass_production_dt));
  }
  // The equation the two above share.
  S establishment_probability(const TF24_Environment<S>& environment,
                              const S& net_mass_production_dt_);

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  S compute_competition(const S& z, const S& height) const;
  // Optimised overload called from Individual<TF24>::compute_competition with the
  // cached competition_effect (= area_leaf(height)) and height_inverse (= 1/height)
  // aux values, matching the shared individual.h interface (no recompute per call).
  S compute_competition(const S& z, const S& area_leaf_,
                        const S& height_inverse) const;
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  S compute_competition(const S& z, const Internals<S>& vars) const {
    return compute_competition(z, vars.aux(aux_idx_competition_effect),
                               vars.aux(aux_idx_height_inverse));
  }

  // The competition contribution and its vertical derivative from one pass, so
  // u^eta is evaluated once. `value` is bit-for-bit the one
  // compute_competition() returns: both read the shading model's own profile,
  // which is what keeps them equal under a flat-top one.
  with_slope<S> compute_competition_and_slope(const S& z, const Internals<S>& vars) const {
    const S& area_leaf_ = vars.aux(aux_idx_competition_effect);
    const S height_inverse = vars.aux(aux_idx_height_inverse);
    const S scale = pars.k_I * area_leaf_;
    const std::pair<S, S> Qq =
      canopy_shape.Q_and_q(z * height_inverse, z, height_inverse);
    // Negated on the way out: q is -dQ/dz and everything above here carries the
    // signed slope.
    return {scale * Qq.first, -(scale * Qq.second)};
  }


  // The fraction of root mass below soil depth `z`, for a plant rooted to
  // `rooting_depth` with shape exponent `eta_x` (pars.root_depth_shape_eta). The
  // canopy's own cumulative form is CanopyShape::Q, at pars.eta.
  S Q(const S& z, const S& rooting_depth, const S& eta_x) const;

  // The inverse of dheight_darea_leaf, so the allometry has one source.
  S darea_leaf_dheight(const S& area_leaf) const {
    return 1.0 / dheight_darea_leaf(area_leaf);
  }

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // The seed's height and leaf area at the current scalar.
  //
  // Preparation solves the height in plain arithmetic and cannot run at an active
  // scalar, so on a differentiated path the height is declared by the residual that
  // defines it -- live mass equals seed mass -- and the leaf area is derived from
  // the height it returns. The two come from one call because the leaf area's own
  // partials in the allometric constants and its chain through the height are the
  // same channel: taking either with the other fixed mixes them.
  struct SeedGeometry { S height; S area_leaf; };
  SeedGeometry seed_geometry() const {
    if constexpr (std::is_same_v<S, double>) {
      return {height_0, area_leaf_0};
    } else {
      // ⚠️ THE SLOPE IS READ, NOT TAKEN HERE. Taking it needs the strategy at a
      // tangent, and a rebind CONSTRUCTS one -- whose Leaf member builds two
      // vulnerability curves from an incomplete gamma per knot, which
      // assign_from then overwrites with a copy. This runs per rate evaluation,
      // so those curves were built and discarded about a hundred thousand times
      // a gradient. prepare_strategy takes the tangent once instead.
      const S h = odelia::implicit_value<S>(
        height_0, dmass_dheight_0,
        [this](const S& y) -> S {
          return mass_live_given_height(y) - pars.omega;
        });
      return {h, area_leaf(h)};
    }
  }

  // Set constants within TF24_Strategy
  void prepare_strategy();

  // The same strategy at scalar U.

  // Another strategy's values, written into this one. prepare_strategy() is
  // refused at an active scalar, so what it produced is carried rather than
  // rebuilt: the leaf, and the canopy shape below.
  //
  // A rebind hands back a fresh strategy, so anything it leaves default has to be
  // set here rather than left alone -- an assignment writes into one that already
  // exists, and a member it does not write keeps a previous state's working.
  template <class S1>
  void assign_from(const TF24_Strategy<S1>& src) {
    // Qualified: these are the base's, and an unqualified name is not looked up
    // in a dependent base.
    this->birth_rate_x = src.birth_rate_x;
    this->birth_rate_y = src.birth_rate_y;
    this->is_variable_birth_rate = src.is_variable_birth_rate;
    this->collect_all_auxiliary = src.collect_all_auxiliary;
    this->size_0 = src.size_0;
    this->control = src.control;
    this->name = src.name;
    this->extrinsic_drivers = src.extrinsic_drivers;

    TF24_Pars<S1> from_pars = src.pars;
    std::vector<S1*> from = from_pars.field_ptrs();
    std::vector<S*> to = pars.field_ptrs();
    // Both lists come from the same table at two scalars, so this says the
    // rebind still describes the same struct.
    util::check_length(from.size(), to.size());
    util::check_length(to.size(), TF24_Pars<S>::field_count);
    for (size_t i = 0; i < from.size(); ++i) {
      *to[i] = S(odelia::util::to_passive(*from[i]));
    }

    shading_model_ = src.shading_model_;
    eta_c = S(odelia::util::to_passive(src.eta_c));
    birth_rate_ = S(odelia::util::to_passive(src.birth_rate_));
    canopy_shape.initialise(pars.eta, shading_model_);
    height_0 = src.height_0;
    dmass_dheight_0 = src.dmass_dheight_0;
    area_leaf_0 = S(odelia::util::to_passive(src.area_leaf_0));
    leaf = src.leaf;
    storage_gate_width = src.storage_gate_width;
    storage_prod_eps = src.storage_prod_eps;
    plasticity_gate_width = src.plasticity_gate_width;
    beta_R_H = src.beta_R_H;
    beta_R_V = src.beta_R_V;
    function_integrator = src.function_integrator;
    crown = src.crown;
    // Shared, not copied: a clamp on this path fires inside a per-unit copy that
    // is then discarded, so the count has to land in storage the run still owns.
    // The missing-row flag is shared for the same reason.
    clamps.differentiated = src.clamps.differentiated;
    curvature_margin = src.curvature_margin;
    recorded_refusal = src.recorded_refusal;
    leaf_points = src.leaf_points;

    // Sized, not copied: these hold one right-hand side's working, and are
    // rewritten before they are read.
    root_carbon_per_leaf_area_.assign(src.root_carbon_per_leaf_area_.size(),
                                      S(0.0));

    // The index maps and the slot numbers are a function of the names, so they
    // are derived rather than carried.
    refresh_indices();
  }

  // This strategy at another scalar, carrying no derivative -- assign_from takes
  // every parameter through to_passive. It is what an implicit node asks a
  // parameter set for when it needs a slope with the parameters fixed.
  template <class U>
  TF24_Strategy<U> rebind() const {
    return TF24_Strategy<U>(tf24_strategy_from_t{}, *this);
  }

  // Birth height of a (germinated) seed. Strategy-agnostic accessor used by
  // the templated Individual; here height_0 is derived in prepare_strategy().
  double initial_height() const { return height_0; }

  // Crown shading model, resolved once from control.shading_model in
  // prepare_strategy(). TF24 supports deep-crown, mean-light (its default)
  // and crown-centre; PPA is not available for TF24.
  ShadingModel shading_model_ = ShadingModel::MeanLight;

  // Biological (user-settable) parameters; see TF24_Pars above.
  TF24_Pars<S> pars;

  // The exponent of the stem path integral, beta = 2*D_c + theta_c. The factor 2
  // on D_c is the packing limit: under a conserved lumen fraction, widening is
  // paid for by proportionally fewer conduits, so sapwood-specific conductivity
  // scales as D^2 and not the D^4 of Hagen-Poiseuille. See
  // plant/stem_hydraulics.h.
  //
  // What it implies: leaf-specific resistance grows as H^(1-beta) rather than
  // linearly with height. beta = 0 is the linear case; the default beta = 0.4
  // gives H^0.6; beta >= 1 saturates, so resistance approaches a finite limit no
  // matter how tall the plant grows. Larger beta therefore means a weaker height
  // penalty on carbon gain. Named for the path rather than for either parameter,
  // since it belongs to neither.
  //
  // DO NOT cache this in a member. rebind() runs before the directions are
  // set on ad_parameters(), so a copy taken at rebind carries the passive
  // value and D_c's column comes back an exact zero -- a number that reads as an
  // answer. Read here, it carries whatever tape identity pars.D_c has at the
  // call. prepare_strategy() validates (beta, L_tip) off this same expression,
  // so nothing reaches stem_hydraulics::effective_path_length that the guards
  // there have not seen.
  S stem_path_exponent() const { return S(2.0) * pars.D_c + pars.theta_c; }

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  S eta_c     = NA_REAL; // crown shape factor, precomputed from pars.eta
  CanopyShape<S> canopy_shape;
  // Height and leaf area of a (germinated) seed
  double height_0  = NA_REAL;
  // dM/dh at the seed height above, which is what closes the implicit function
  // theorem on it in seed_geometry. Derived here rather than there because the
  // tangent that takes it needs a whole strategy rebound, and seed_geometry runs
  // on the rate-evaluation path. NA_REAL until prepare_strategy runs, so a
  // strategy that skipped it stops in implicit_value's own dF/dy guard rather
  // than returning a number.
  double dmass_dheight_0 = NA_REAL;
  S area_leaf_0;

  // Embedded leaf hydraulic/photosynthesis sub-model, built in prepare_strategy()
  Leaf leaf;

  // Width of the smooth reserve gate G(r) on growth (#517); small relative to
  // [0,1] so the switch about the growth threshold a_st2 is fairly sharp but
  // differentiable. storage_prod_eps smooths the positive-part of net production
  // (replacing the old hard net>0 growth cutoff) for AD-readiness.
  double storage_gate_width = 0.1;
  double storage_prod_eps   = 1e-4;
  // Width of the replacement gate, in units of the marginal leaf return (#516).
  double plasticity_gate_width = 0.02;

  // The clamp sites this strategy reaches; the list itself is shared with the
  // environment, which reaches the others.
  clamp_counter clamps;

  // Count a clamp against the path it fired on. Which path this is is a property
  // of the scalar, so it is decided at compile time.
  void note_clamp(int site) {
    if constexpr (std::is_same_v<S, double>) {
      ++clamps.forward[site];
    } else {
      ++(*clamps.differentiated)[site];
    }
  }

  // How close the interior derivation came to dividing by nothing. The guard
  // below refuses on a floor, and a floor nothing approaches and a floor nothing
  // reaches report the same green -- so the margin is carried out beside the
  // count rather than left to be assumed. Shared for the same reason the counts
  // are: the copy that measures it is discarded.
  std::shared_ptr<double> curvature_margin = std::make_shared<double>(-1.0);
  void note_curvature(double value) {
    const double m = std::abs(value);
    if (*curvature_margin < 0.0 || m < *curvature_margin) {
      *curvature_margin = m;
    }
  }

  // Below this the collar's own response is amplification rather than an answer.
  // It lives on Control because it moves which states answer, so two gradients
  // taken at different values are gradients of different functions and
  // stand_gradient() has to be able to say so.
  //
  // The units are the profit's, so a reparameterisation that rescales profit
  // needs it re-measured -- which is what curvature_margin is for.
  double curvature_floor() const { return this->control.gradient_curvature_floor; }

  // Where a row a census metric needs does not exist. The values still do, so the
  // recording carries them on and leaves the missing rows off the tape.
  //
  // Recorded rather than thrown, and a recording is the reason: a throw from here
  // unwinds through a live tape and the whole sweep to deliver what a poll after
  // it delivers, and there is nothing for it to save, because a metric is a sum
  // and a sum has no defined value with an undefined term.
  //
  // Shared for the reason the clamp counts are: the strategy that records it is a
  // per-unit copy the sweep discards, so a refusal on the copy is one nothing can
  // read. Cleared where the call that reads it starts, so it latches for one
  // gradient call rather than for one operating point.
  std::shared_ptr<refusal> recorded_refusal = std::make_shared<refusal>();

  // How many leaf solves landed in each operating-point kind, indexed by the
  // enum. Diagnostic and deliberately outside rebind: a block copies the
  // strategy per unit and discards it, so carrying the tally across would count
  // the sweep's copies as well as the run.
  std::vector<size_t> operating_point_counts =
    std::vector<size_t>(Leaf::operating_point_kind_count, 0);

  // The root-architecture model's two constants. The leaf takes resistances and
  // this strategy owns the model that produces them: the 1/3 : 2/3 root split
  // and the dz^2 vertical scaling are not gas exchange.
  double beta_R_H = 3.4e2;
  double beta_R_V = 9.4e3;

  // Cached aux/state indices, resolved once in refresh_indices(), so the hot
  // compute_rates path does not do a std::map<string,int>::at (string compare)
  // lookup per ODE derivs evaluation per individual (profile hot spot).
  int aux_idx_competition_effect = -1;
  int aux_idx_height_inverse = -1;
  int aux_idx_net_mass_production_dt = -1;
  int aux_idx_root_mass = -1;
  int aux_idx_opt_psi_stem = -1;
  int aux_idx_opt_root_psi = -1;
  int aux_idx_transpiration = -1;
  int aux_idx_E_up = -1;
  int aux_idx_profit = -1;
  int aux_idx_shadow_cost = -1;
  int aux_idx_stom_cond_CO2 = -1;
  int aux_idx_assimilation = -1;
  int aux_idx_Tleaf = -1;
  int aux_idx_leaf_marginal_return = -1;
  int aux_idx_sapwood_marginal_return = -1;
  int aux_idx_area_sapwood = -1;       // only present when collect_all_auxiliary
  int state_idx_area_heartwood = -1;
  int state_idx_mass_heartwood = -1;
  int state_idx_storage        = -1;
  // ln(area_leaf / area_leaf(height)): how far the canopy sits from the leaf
  // area the height allometry prefers. Signed, so not a non_negative_state.
  int state_idx_log_area_leaf_departure = -1;
  // ln(area_sapwood / (theta * area_leaf)): how far the conducting stem sits
  // from the leaf area it supports, i.e. the departure of the Huber value from
  // its preferred one. Also signed.
  int state_idx_log_area_sapwood_departure = -1;
  // Running mean of the marginal leaf's carbon balance. The gate reads THIS,
  // not the instantaneous value, so a plant does not shed and re-flush with
  // every passing dry spell.
  int state_idx_leaf_marginal_tracked = -1;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;
  // The quadrature over this strategy's unit crown, which the mean light reads
  // the light field on.
  crown_rule crown;

  // Reusable per-layer root-carbon buffer, refilled (not reallocated) each
  // net_mass_production_dt call to avoid a heap allocation per derivs eval.
  // Carries S: per-layer root carbon is one of the state directions the leaf's
  // supplied Jacobian has rows for, so it is a live gradient channel.
  std::vector<S> root_carbon_per_leaf_area_;

  // The same numbers with the derivative stripped, because the architecture
  // model and the leaf both take double. Held rather than made per call for the
  // reason the buffer above is.
  std::vector<double> root_carbon_value_;

  // Filled in place by the architecture model each call and taken by the leaf as
  // a const reference, so it must not be moved from: that is what keeps its
  // capacity across calls.
  phylloptim::RootNetwork root_network_;

  // The soil potentials with the derivative stripped, for the reason
  // root_carbon_value_ is.
  std::vector<double> psi_soil_value_;

  // The first reason is kept: a later one is a consequence of the same
  // degeneracy.
  void refuse(const std::string& why) {
    if (!recorded_refusal->happened()) {
      recorded_refusal->reason = why;
    }
  }

  // The leaf's own sites, one site at a time. Reading them as a folded vector
  // allocates, and a differentiated leaf evaluation reads them twice.
  using leaf_clamp_tally =
      std::array<std::size_t, phylloptim::CLAMP_SITE_COUNT>;
  static_assert(CLAMP_LEAF_FIRST + phylloptim::CLAMP_SITE_COUNT ==
                    CLAMP_SITE_COUNT,
                "the leaf's sites are the last of plant's, so this tally maps "
                "onto them by adding CLAMP_LEAF_FIRST");

  leaf_clamp_tally leaf_clamps() const {
    leaf_clamp_tally out{};
    for (std::size_t s = 0; s < out.size(); ++s) {
      out[s] = leaf.clamp_count(static_cast<int>(s));
    }
    return out;
  }

  // What the leaf clamped between two readings of its own tally. Everything
  // between re-supplies and re-solves it many times, so anything clamped there
  // is the gradient's and the forward run's share is the total less this.
  void note_leaf_clamps(const leaf_clamp_tally& before) {
    const leaf_clamp_tally after = leaf_clamps();
    for (std::size_t s = 0; s < after.size(); ++s) {
      if (after[s] > before[s]) {
        (*clamps.differentiated)[CLAMP_LEAF_FIRST + s] += after[s] - before[s];
      }
    }
  }

  // What this species' leaf solves found, kept against the rate evaluation that
  // found them so a pass re-running the model over the same states places the point
  // instead of searching for it again. The patch hands the address down.
  //
  // Shared, not copied, for the reason the clamp counts are: the pass that reads the
  // record runs on a rebound strategy.
  std::shared_ptr<leaf_solved_points> leaf_points =
      std::make_shared<leaf_solved_points>();
  // What one of this species' rate evaluations solves for.
  using solved_values = leaf_solved_points::points;
  void store_solved(solved_values& into) { leaf_points->store_into(into); }
  void load_solved(const solved_values& from) { leaf_points->load_from(from); }
  void end_solved() { leaf_points->end(); }
  size_t leaf_placements() const { return leaf_points->placements(); }


  // The leaf's two outputs on the active chain, carrying its supplied Jacobian.
  // Written by net_mass_production_dt before compute_rates reads either.
  S leaf_profit_;
  std::vector<S> leaf_soil_consumption_;
  // The shadow price of the water the leaf used, `lambda_o * E`: added back to
  // the objective to give the carbon kept (net_mass_production_dt says why).
  S leaf_shadow_cost_;

  // Flexible allometry's working (#516), written by net_mass_production_dt and
  // read by compute_rates in the same evaluation.
  //
  // The marginal leaf's carbon balance; dimensionless, zero at break-even.
  S marginal_leaf_return_ = 0.0;
  // The withheld shares the budget actually charged, kept so compute_rates
  // cannot form them differently from the budget that used them. The carbon not
  // spent and the tissue lost are two halves of one decision, so they must come
  // from one number.
  S withheld_leaf_ = 0.0;
  S withheld_sapwood_ = 0.0;
  // Respiration plus turnover (kg yr^-1), the strictly-positive scale the
  // sapwood controller normalises by.
  S maintenance_flux_ = 0.0;
  // The leaf-specific max hydraulic conductance actually handed to the leaf
  // (kg m^-2 LA s^-1 MPa^-1), stored so the sapwood controller differentiates
  // the value the model used rather than a recomputed copy of it.
  S leaf_specific_conductance_max_ = 0.0;
  // d(profit)/d(kmax) at the leaf's own operating point, umol m^-2 s^-1 per unit
  // conductance. A finite difference through the leaf, so a double; zero unless
  // sapwood acclimation is on.
  double dprofit_dkmax_ = 0.0;

  // Every active value this strategy holds: the whole parameter table, the
  // constant birth rate, the leaf outputs, and the per-layer root carbon. `leaf`
  // is not among them -- it solves in double, which is what keeps the leaf off
  // the tape.
  // ⚠️ A LEAF OUTPUT LEFT OFF THIS LIST CARRIES NO ROWS, SILENTLY.
  template <class F>
  void for_each_active(F&& f) {
    for (const ad_parameter& field : pars.ad_parameter_table()) {
      f(pars.*field.at);
    }
    f(birth_rate_);
    odelia::ode::visit_active(f, leaf_profit_, leaf_soil_consumption_,
                              leaf_shadow_cost_, root_carbon_per_leaf_area_,
                              marginal_leaf_return_, withheld_leaf_,
                              withheld_sapwood_, maintenance_flux_,
                              leaf_specific_conductance_max_);
  }
};

template <typename S>
typename TF24_Strategy<S>::ptr make_strategy_ptr(TF24_Strategy<S> s);

// Named here for clarity; promotion to user-tunable traits (RcppR6) is a
// deliberate follow-up (see vignettes/models/code_review_leaf_tf24.qmd #9).
// rescales total fine-root mass into the per-layer carbon units expected by the
// root hydraulic network in Leaf::set_physiology.
inline const double root_mass_carbon_scale = 83.26 * 0.5;
// ⚠️ DO NOT HOIST the per-second -> annual factor 60*60*12*365 (seconds of
// daylight per year, 12 h day x 365 d) to a constant. It recurs in compute_rates
// and net_mass_production_dt below; collapsing the 4-step integer product into
// one double changes the floating-point rounding, and the adaptive ODE amplifies
// it (offspring_production shifts ~0.2%).

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
template <typename S>
TF24_Strategy<S>::TF24_Strategy() {
  this->collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  this->name = "TF24";
}

// assumes optimise_psi_stem_TF has been run for optimal psi_stem
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
                                           const S& conductance_max) {
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
  // Transpiration at the collar the solve placed, carrying its move: the draw is
  // taken at the passive collar, so its flux alone would miss dcollar/dtheta. The
  // step is zero in value and its slope is the draw's own. Read for the shadow
  // cost below.
  S transpiration_at = S(leaf.transpiration_);
  try {
    if (leaf.operating_point_kind() == Leaf::OperatingPointKind::Interior) {
      slope = leaf.template profit_curvature_at_collar<
          phylloptim::Leaf::CostCurve::TF24_floor>();
      // The interior derivation divides by the profit's curvature. At a bound
      // the slope is the bound's own and this floor is not about it.
      note_curvature(slope);

      // ⚠️ THE DISAGREEMENT PROBE THAT SAT HERE IS GONE, AND SO IS WHAT IT
      // CHECKED. It differenced the marginal and compared that against the
      // curvature, because the curvature was an ANALYTIC assembly that could be
      // wrong in a way nothing else could see -- on the century stand it caught a
      // +34.41 where the difference said -9.63.
      //
      // profit_curvature_at_collar IS that difference now, measured rather than
      // assembled, so the probe would compare a difference with itself. The
      // reason it is a difference is in phylloptim: the analytic route needs four
      // significant digits held through a 1e+05 cancellation and then multiplied
      // by 8.1e+04, and returns an eighth of the answer with the right sign.
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
    transpiration_at =
        draw.E_up.value + draw.E_up.slope * (collar - S(leaf.opt_root_psi_));
  } catch (const std::runtime_error& e) {
    note_leaf_clamps(clamps_before);
    // The leaf produced nothing to record a row against, so every output takes
    // the value its own solve fixed and carries none -- which is what a refusal
    // costs everywhere else here, and lets the sweep reach the poll that reads
    // this.
    refuse(std::string("TF24 gradient: ") + e.what());
    leaf_profit_ = S(leaf.profit_);
    leaf_shadow_cost_ = S(leaf.shadow_cost());
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

  carry(leaf.profit_, got.profit, leaf_profit_, true);
  leaf_soil_consumption_.assign(n_layer, S(0.0));
  util::check_length(got.uptake.size(), n_layer);
  for (std::size_t j = 0; j < n_layer; ++j) {
    carry(leaf.soil_consumption_[j], got.uptake[j], leaf_soil_consumption_[j],
          false);
  }
  // The shadow cost, recorded against the leaf's own value like every output
  // here. It is `lambda_o * E` on TF24_floor and zero on any other curve -- the
  // same condition phylloptim's shadow_cost() applies -- so a curve that reads no
  // price carries no row through it. It reads the water, so a refused water row
  // leaves it without one too.
  const S shadow_rows =
      leaf.cost_curve_ == Leaf::CostCurve::TF24_floor
          ? in[phylloptim::par_TF24_floor_lambda_o] * transpiration_at
          : S(0.0);
  carry(leaf.shadow_cost(), shadow_rows, leaf_shadow_cost_, false);
}

template <typename S>
void TF24_Strategy<S>::refresh_indices () {
    // Create and fill the name to state index maps
  this->state_index = std::map<std::string,int>();
  this->aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (size_t i = 0; i < state_names_vec.size(); i++) {
    this->state_index[state_names_vec[i]] = i;
  }
  for (size_t i = 0; i < aux_names_vec.size(); i++) {
    this->aux_index[aux_names_vec[i]] = i;
  }

  // Cache integer indices for the keys used in the hot compute_rates path, so it
  // does no std::map<string,int> lookup per derivs evaluation.
  aux_idx_competition_effect    = this->aux_index.at("competition_effect");
  aux_idx_height_inverse        = this->aux_index.at("height_inverse");
  aux_idx_net_mass_production_dt = this->aux_index.at("net_mass_production_dt");
  aux_idx_root_mass             = this->aux_index.at("root_mass");
  aux_idx_opt_psi_stem          = this->aux_index.at("opt_psi_stem");
  aux_idx_opt_root_psi          = this->aux_index.at("opt_root_psi");
  aux_idx_transpiration         = this->aux_index.at("transpiration");
  aux_idx_E_up                  = this->aux_index.at("E_up_");
  aux_idx_profit                = this->aux_index.at("profit");
  aux_idx_shadow_cost           = this->aux_index.at("shadow_cost");
  aux_idx_stom_cond_CO2         = this->aux_index.at("stom_cond_CO2");
  aux_idx_assimilation          = this->aux_index.at("assimilation");
  aux_idx_Tleaf                 = this->aux_index.at("Tleaf");
  aux_idx_leaf_marginal_return  = this->aux_index.at("leaf_marginal_return");
  aux_idx_sapwood_marginal_return =
    this->aux_index.at("sapwood_marginal_return");
  // area_sapwood is only registered when collect_all_auxiliary is set.
  aux_idx_area_sapwood = this->aux_index.count("area_sapwood") ? this->aux_index.at("area_sapwood") : -1;
  state_idx_area_heartwood      = this->state_index.at("area_heartwood");
  state_idx_mass_heartwood      = this->state_index.at("mass_heartwood");
  state_idx_storage             = this->state_index.at("storage");
  state_idx_log_area_leaf_departure =
    this->state_index.at("log_area_leaf_departure");
  state_idx_log_area_sapwood_departure =
    this->state_index.at("log_area_sapwood_departure");
  state_idx_leaf_marginal_tracked =
    this->state_index.at("leaf_marginal_return_tracked");
  check_state_layout(this->state_index, state_size(), "TF24");
}

// [eqn 2] area_leaf (inverse of [eqn 3])
template <typename S>
S TF24_Strategy<S>::area_leaf(const S& height) const {
  return pow(height / pars.a_l1, 1.0 / pars.a_l2);
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
template <typename S>
S TF24_Strategy<S>::mass_leaf(const S& area_leaf) const {
  return area_leaf * pars.lma;
}

// [eqn 4] area and mass of sapwood
template <typename S>
S TF24_Strategy<S>::area_sapwood(const S& area_leaf) const {
  return area_leaf * pars.theta;
}

// The conducting area a plant actually has, as against the one its leaf area
// prefers. At a resting departure sapwood_per_leaf_area() returns pars.theta
// unchanged, so this is bit-identical to the overload above.
template <typename S>
S TF24_Strategy<S>::area_sapwood(const S& area_leaf,
                                 const S& sapwood_departure) const {
  return area_leaf * sapwood_per_leaf_area(sapwood_departure);
}

template <typename S>
S TF24_Strategy<S>::mass_sapwood(const S& area_sapwood, const S& height) const {
  return area_sapwood * height * eta_c * pars.rho;
}

// [eqn 5] area and mass of bark
template <typename S>
S TF24_Strategy<S>::area_bark(const S& area_leaf) const {
  return pars.a_b1 * area_leaf * pars.theta;
}

template <typename S>
S TF24_Strategy<S>::mass_bark(const S& area_bark, const S& height) const {
  return area_bark * height * eta_c * pars.rho;
}

template <typename S>
S TF24_Strategy<S>::area_stem(const S& area_bark, const S& area_sapwood,
                            const S& area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

template <typename S>
S TF24_Strategy<S>::diameter_stem(const S& area_stem) const {
  using std::sqrt;
  return sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
template <typename S>
S TF24_Strategy<S>::mass_root(const S& area_leaf) const {
  return pars.a_r1 * area_leaf;
}

// [eqn 7b] Mass of coarse (structural) roots
template <typename S>
S TF24_Strategy<S>::mass_coarse_root(const S& mass_sapwood) const {
  return pars.a_cr1 * mass_sapwood;
}

// [eqn 8] Total mass
// The coarse-root term is appended, not interleaved, so that with the default
// a_cr1 = 0 the sum is bit-identical to the pre-#349 expression (x + 0.0 == x
// exactly, and the association order of the other terms is untouched).
template <typename S>
S TF24_Strategy<S>::mass_live(const S& mass_leaf, const S& mass_bark,
                           const S& mass_sapwood, const S& mass_root,
                           const S& mass_coarse_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root + mass_coarse_root;
}

template <typename S>
S TF24_Strategy<S>::mass_total(const S& mass_leaf, const S& mass_bark,
                            const S& mass_sapwood, const S& mass_heartwood,
                            const S& mass_root,
                            const S& mass_coarse_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root +
    mass_coarse_root;
}

template <typename S>
S TF24_Strategy<S>::mass_above_ground(const S& mass_leaf, const S& mass_bark,
                            const S& mass_sapwood, const S& mass_heartwood) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
}

// for updating auxiliary state
template <typename S>
void TF24_Strategy<S>::update_dependent_aux(const int index, Internals<S>& vars) {
  // Leaf area depends on BOTH height and the departure state, so either write
  // has to refresh it. Individual::set_ode_state refreshes per slot in slot
  // order, so the height pass (slot 0) necessarily reads the departure still
  // sitting in vars from the previous evaluation; the departure pass recomputes
  // from both, which is why it must not be a partial update. Whichever fires
  // last leaves the aux consistent with the whole state vector.
  if (index == HEIGHT_INDEX || index == state_idx_log_area_leaf_departure) {
    const S& height = vars.state(HEIGHT_INDEX);
    // A = A*(h) * exp(phi). With plasticity off phi is exactly 0 and exp(0) is
    // exactly 1, so this is bit-identical to the plain area_leaf(height) it
    // replaces -- the multiplication is exact, and there is no addition for an
    // FMA to fuse.
    vars.set_aux(aux_idx_competition_effect, area_leaf_carried(vars));
    vars.set_aux(aux_idx_height_inverse, 1.0 / height);
  }
}


// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
template <typename S>
void TF24_Strategy<S>::compute_rates(const TF24_Environment<S>& environment,  Internals<S>& vars) {
  const S& height = vars.state(HEIGHT_INDEX);
  const S& area_leaf_ = vars.aux(aux_idx_competition_effect);
  // The Huber value's departure from theta. Read once here and threaded through
  // every sapwood quantity below, so the conducting area a plant has and the
  // one its leaf area prefers cannot disagree between two call sites.
  const S& sapwood_departure = vars.state(state_idx_log_area_sapwood_departure);
  const S& leaf_departure = vars.state(state_idx_log_area_leaf_departure);

  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_,
                           vars.aux(aux_idx_height_inverse),
                           sapwood_departure, leaf_departure,
                           vars.state(state_idx_leaf_marginal_tracked));

  // The running mean chases the instantaneous balance at a_pl4. A first-order
  // filter, so the memory is 1/a_pl4 and nothing about it can overshoot.
  vars.set_aux(aux_idx_leaf_marginal_return, marginal_leaf_return_);
  vars.set_rate(state_idx_leaf_marginal_tracked,
                pars.a_pl4 * (marginal_leaf_return_ -
                              vars.state(state_idx_leaf_marginal_tracked)));

  // Read back the shares the budget above actually charged, rather than forming
  // them again: the carbon not spent and the tissue lost are two halves of one
  // decision, and they balance only if they come from one number.
  const S withheld         = withheld_leaf_;
  const S withheld_sapwood = withheld_sapwood_;
  const S replacement      = 1.0 - withheld;

  // Read by establishment_probability for the boundary node, so this one is a
  // rate's input rather than a reading of it.
  vars.set_aux(aux_idx_net_mass_production_dt, net_mass_production_dt_);

  // The rest are readings: their only consumers are r_internals and ode_aux,
  // both R-facing. Writing them at double alone was tried and measured at no
  // difference, because a store into a slot the vector already holds pushes a
  // statement and registers nothing -- it is a newly constructed active scalar
  // that costs a slot, not a write to a live one.
  vars.set_aux(aux_idx_root_mass, mass_root(area_leaf_));
  vars.set_aux(aux_idx_opt_psi_stem, leaf.opt_psi_stem_);
  vars.set_aux(aux_idx_opt_root_psi, leaf.opt_root_psi_);
  vars.set_aux(aux_idx_transpiration, leaf.transpiration_);
  vars.set_aux(aux_idx_E_up, leaf.E_up_);
  vars.set_aux(aux_idx_profit, leaf.profit_);
  // ⚠️ THE INDEX WAS DECLARED AND NEVER ASSIGNED, AND THE AUX NEVER WRITTEN, so
  // shadow_cost reported whatever the slot last held -- zero at the default
  // price, which is also its correct value there, so nothing could tell the two
  // apart until a price was set.
  vars.set_aux(aux_idx_shadow_cost, S(leaf.shadow_cost()));
  vars.set_aux(aux_idx_stom_cond_CO2, leaf.stom_cond_CO2_);
  vars.set_aux(aux_idx_assimilation, leaf.assim_colimited_);
  // The leaf's own temperature at the operating point, not the `leaf_temp`
  // driver -- see aux_names(). Equal to the driver while pars.use_energy_balance
  // is off; solved from the leaf's transpiration when it is on.
  vars.set_aux(aux_idx_Tleaf, leaf.Tleaf_);




  // consumption rates should be emerging from net_mass_produciton_dt
  // convert evapotranspiration per leaf area per soil layer (mol H20 m^-2 s^-1) to canopy-level total 
  // yearly evapotranspiration per soil layer (m yr^-1)
  // stubbing out E_p for integration
  int soil_number_of_depths_ = environment.get_soil_number_of_depths();


  for (int i = 0; i < soil_number_of_depths_; i++) {

    // evapotranspiration (mol H20 m^-2 s^-1 layer^-1)
    // consumption rate (m yr^-1 layer ^-1)
    vars.set_consumption_rate(i, evapotranspiration_dt(area_leaf_, i)*60*60*12*365/1000*kg_per_mol_h2o);
  }

  // A storage pool buffers demography against short-term productivity swings:
  // growth is gated on having ample reserves, and mortality reads the buffered
  // relative reserves rather than instantaneous net production -- so a trough
  // draws reserves down and death is gradual, instead of the growth cutoff and
  // ~1e32 mortality spike that caused the #550 blow-up.
  // Read unclamped: a stage below empty keeps finite rates for the error
  // estimate, and non_negative_states() refuses a step that ends there.
  const S storage     = vars.state(state_idx_storage);
  const S storage_max = storage_capacity(area_leaf_, height, sapwood_departure);
  using std::exp;
  using std::sqrt;
  const S r           = storage_max > 0.0 ? storage / storage_max : S(0.0);
  // Reserve-gated growth (#517), following Daniel's intuition that a plant
  // should not grow unless it has ample carbon in storage. Growth and
  // reproduction proceed at the *production* rate, but scaled by a smooth gate
  // G(r) that is ~0 at low relative reserves and ~1 near capacity (logistic
  // centred on the growth threshold a_st2). Carbon not spent on growth refills
  // storage (dS/dt below), so: full reserves -> grow at production rate (healthy
  // dynamics preserved); low reserves -> redirect carbon to refill, growth
  // pauses; drought (net<0) -> reserves drain, growth halts and only resumes
  // once they refill (buffered growth AND mortality). Note growth is gated at
  // the production rate rather than metered out of the pool as a_st2*S: with
  // sapwood-scaled capacity, storage is tiny relative to seedling productivity,
  // so a discharge-rate limit would choke establishment.
  const S G = 1.0 / (1.0 + exp(-(r - pars.a_st2) / storage_gate_width));
  // Smooth positive part of net production (replaces the old hard net>0 cutoff).
  const S P = net_mass_production_dt_;
  const S Ppos =
    0.5 * (P + sqrt(P * P + storage_prod_eps * storage_prod_eps));
  const S growth_flux = Ppos * G;

  const S fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
  const S fraction_allocation_growth_ = fraction_allocation_growth(height);

  // --- Where the growth flux goes (#516) ------------------------------------
  // Growth moves a plant ALONG its allometry, plasticity moves it OFF, and
  // rebuilding moves it back. Keeping those three separate is what makes the
  // resting model exact and the departed one interpretable.
  //
  // So the allometric derivatives are evaluated at the leaf area the plant's
  // HEIGHT implies, area_leaf(height), not at the leaf area it actually has.
  // They are properties of the curve, and a plant's place on the curve is given
  // by its height; its leaf area is a departure from that. Evaluated this way,
  // extension is exactly departure-preserving for any departure, so the whole
  // of the departure's motion is plasticity. At rest the two coincide bitwise
  // (exp(0) == 1), so this is the same number the fixed model used.
  const S area_leaf_star = area_leaf(height);
  const S darea_leaf_dmass_live_ =
    darea_leaf_dmass_live(area_leaf_star, sapwood_departure);

  // The share of growth redirected into closing the canopy gap: it vanishes
  // exactly when the gap closes, widens as the gap does, and is gated by the
  // SAME signal that decides replacement.
  //
  // ⚠️ THE `replacement` FACTOR IS LOAD-BEARING, not a refinement. Without it
  // the gap alone drives rebuilding, so a starving plant withholds replacement
  // and spends its growth flux rebuilding at the same time -- contradictory,
  // since carbon enough to rebuild is carbon enough to maintain. Measured, the
  // unconditioned form rebuilt at ~1.5/yr against shedding's 0.42/yr, so the
  // canopy thinned by 5 per cent and stalled, and the Huber value fell instead
  // of rising because rebuilding adds leaf area without sapwood. Gating both on
  // the same signal makes a plant rebuild exactly to the degree it is willing
  // to maintain.
  //
  // The departure cannot go positive under these dynamics -- shedding only
  // lowers it and this share only raises it, stopping at zero -- so the form
  // needs no bound.
  const S rebuild_share =
    replacement * (1.0 - exp(leaf_departure / pars.a_pl2));
  const S extension_share = 1.0 - rebuild_share;

  const S area_leaf_dt =
    extension_share * growth_flux * fraction_allocation_growth_ *
      darea_leaf_dmass_live_;

  // Height follows extension alone, so rebuilding a canopy buys no height and
  // height still cannot fall: every factor here is non-negative.
  vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_star) * area_leaf_dt);
  vars.set_rate(FECUNDITY_INDEX,
    fecundity_dt(growth_flux, fraction_allocation_reproduction_));

  // Rebuilding buys leaf area at unchanged height, and it buys the WHOLE
  // package -- leaf, fine root, and the sapwood and bark that supply it -- at
  // the ratio the plant prefers rather than the one it has. Buying leaf alone
  // diluted the Huber value, and measured, that dilution dominated: plants
  // finished under-built rather than over-built (#516).
  const S rebuild_rel =
    rebuild_share * growth_flux * fraction_allocation_growth_ /
      (mass_live_per_area_leaf_at_height(height) * area_leaf_);

  // Sapwood -> heartwood conversion is turnover-driven, so it proceeds
  // regardless of carbon status. DO NOT gate it behind net production > 0.
  vars.set_rate(state_idx_area_heartwood,
                area_heartwood_dt(area_leaf_, sapwood_departure));
  const S area_sapwood_ = area_sapwood(area_leaf_, sapwood_departure);
  const S mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  vars.set_rate(state_idx_mass_heartwood, mass_heartwood_dt(mass_sapwood_));

  if (this->collect_all_auxiliary) {
    vars.set_aux(aux_idx_area_sapwood, area_sapwood_);
  }

  // Storage dynamics, as a charge and a drain that are each non-negative
  // without a test: sqrt(P^2 + eps^2) >= |P| makes Ppos and the drain both >= 0,
  // and G is a logistic so 1 - G >= 0. Their difference is the net flux exactly,
  // Ppos(1 - G) - (Ppos - P) = P - growth_flux, so the split itself moves
  // nothing; what it buys is somewhere to limit each direction on its own
  // (#609).
  const S charge = Ppos * (1.0 - G);      // surplus the gate withheld
  const S drain  = Ppos - P;              // shortfall met from reserves
  // Each direction is limited by the room the other has: the charge fills
  // headroom, the drain spends contents. So dS/dt = charge - (charge + drain) r,
  // the pool is a first-order filter on production, and both bounds follow from
  // the form with no scale left to choose -- at r = 0 the rate is charge >= 0,
  // at r = 1 it is -drain <= 0. It is also smooth where the net flux changes
  // sign, where the old branch stepped the slope by (r + drain_ref)/r, unbounded
  // as the reserves empty.
  //
  // Narrowing the drain to r/(r + D) is the shape this replaced and the one a
  // reader is most likely to restore, since it looks like the more faithful
  // limiter. It puts an attracting fixed point at r ~ D whose relaxation time is
  // storage_max D / drain -- under an hour for a seedling at D = 1e-3, against a
  // solver stepping in days, and measured at 26 times the accepted steps.
  //
  // The pool is capped by withholding the surplus rather than by spending it,
  // so production and the two flows no longer balance: at capacity the charge
  // the gate withheld, Ppos(1 - G) ~ 1.2e-4 of production, leaves the budget.
  //
  // The flow relaxes the pool over storage_max/(charge + drain), hours in a
  // seedling; the divisor adds storage_relaxation_offset to that time.
  const S relaxation_rate =
      storage_max > 0.0 ? (charge + drain) / storage_max : S(0.0);
  vars.set_rate(state_idx_storage,
                (charge * (1.0 - r) - drain * r) /
                    (1.0 + relaxation_rate * pars.storage_relaxation_offset));

  // --- The two departures (#516) --------------------------------------------
  // Relative rates, so both are scale-free: the same withheld fraction thins a
  // seedling's canopy and an adult's by the same proportion, which is what
  // keeps these parameters from repeating storage_prod_eps's mistake (#620) of
  // being an absolute constant against fluxes spanning six orders.
  //
  // Leaf area loses the turnover it did not replace and regains what rebuilding
  // buys. Thinning is therefore bounded by k_l -- a plant declines to replace
  // what died rather than actively shedding -- so the canopy cannot empty
  // faster than its leaves turn over.
  vars.set_rate(state_idx_log_area_leaf_departure,
                rebuild_rel - withheld * pars.k_l);

  // --- Sapwood chases the ratio that maximises GROWTH -------------------------
  // An integral controller on the marginal return: dpsi/dt = a_sw * R_s, which
  // converges exactly onto R_s = 0 and needs no separate tracked state, because
  // integrating IS averaging.
  //
  // R_s is the marginal return of sapwood area, in growth rather than in
  // production. Growth is F * f_g * D with D = darea_leaf_dmass_live, so
  //     d(growth)/dpsi  proportional to  dP/dpsi - P * (dm_sap/dA) * D
  // the second term being the price of the extra stem, paid in forgone leaf
  // area. Dropping it optimises production instead, which peaks near 3.3x the
  // pipe-model ratio against growth's 1.4-1.8x.
  S sapwood_marginal_ = 0.0;
  if (odelia::util::to_passive(pars.a_sw) > 0.0 &&
      maintenance_flux_ > 0.0) {
    const S mass_sapwood_now = mass_sapwood_;
    // Read, do not recompute: #617 replaced the height-linear denominator with
    // a path integral, and a recomputed copy here would still have been on the
    // old relation -- wrong by a factor 0.35 at 16 m, and silent.
    const S kmax_now = leaf_specific_conductance_max_;
    const S conv = pars.a_bio * pars.a_y;
    // psi enters kmax and sapwood mass both as exp(psi), so d/dpsi of each is
    // itself -- which is what makes these two terms so simple.
    // ⚠️ SAPWOOD TURNOVER IS GATED, so its psi-derivative is not k_s * m_s. What
    // the budget is actually charged is
    //     replacement_sapwood * k_s * m_s,  replacement_sapwood ∝ exp(-psi/a_pl2)
    // against m_s ∝ exp(psi), so the charge falls as exp((1 - 1/a_pl2) * psi) --
    // steeply, a_pl2 being well below one. Differentiating the UNGATED cost
    // instead drops a term worth +6.8 kg/yr of 13.1 at the default, and moved
    // the controller's zero from 1.35x the pipe-model ratio down to 1.05x.
    // Respiration carries no such gate, so it keeps the plain exp(psi).
    const S charged_sapwood_turnover =
      (1.0 - withheld_sapwood) * pars.k_s * mass_sapwood_now;
    const S dP_dpsi =
      conv * area_leaf_ * dprofit_dkmax_ * kmax_now * 60*60*12*365/1e6
      - conv * pars.r_s * mass_sapwood_now
      - charged_sapwood_turnover * (1.0 - 1.0 / pars.a_pl2);
    const S dmass_sap_dA =
      dmass_sapwood_darea_leaf(area_leaf_star, sapwood_departure);
    // ⚠️ SCALED BY CARBON AVAILABILITY, because re-proportioning the stem is an
    // allocation decision and has to be paid for out of the growth flux.
    // Without this the drift raises area_sapwood at fixed leaf area -- real
    // tissue -- while the budget is charged nothing, so a plant in deficit would
    // keep thickening on carbon it does not have. The reserve gate G is NOT the
    // right switch on its own: at soil 0.13 reserves are still full so G ~ 1,
    // and what has gone to zero is Ppos. growth_flux carries both.
    //
    // The factor is a share of throughput, so it is smooth, lies in [0, 1), and
    // vanishes with the growth flux. Symmetric on purpose: scaling only the
    // building direction would put a kink at R_s = 0, and the turnover route
    // out of an over-built stem stays open regardless, being a real mass flow.
    // Being strictly positive it cannot move the zero, so the controller still
    // converges on the growth optimum -- it only sets how fast, which is what
    // a_sw is for.
    //
    // Divided by the maintenance flux to make it dimensionless: strictly
    // positive, unlike P, so the sign of R_s is the sign of the growth
    // derivative everywhere including deep deficit.
    const S availability = growth_flux / (growth_flux + maintenance_flux_);
    sapwood_marginal_ = availability *
      (dP_dpsi - net_mass_production_dt_ * dmass_sap_dA *
         darea_leaf_dmass_live_) / maintenance_flux_;
  }
  vars.set_aux(aux_idx_sapwood_marginal_return, sapwood_marginal_);

  // The Huber value's departure is the DIFFERENCE of the two losses, so the
  // drought acclimation falls out of the arithmetic rather than being imposed:
  // leaf area thins at up to k_l = 0.457/yr while conducting area is only lost
  // to heartwood at k_s = 0.2/yr, and withholding both raises sapwood per leaf
  // area. The excess is then renewed less (withheld_sapwood > withheld), which
  // is the restoring force.
  //
  // The rebuilding term is exact geometry, not a parameter: rebuilding adds leaf
  // area and sapwood in the ratio theta, so it moves the two logs by b/A and
  // theta*b/A_s, whose difference is (b/A) * (exp(-psi) - 1). It vanishes at
  // psi = 0, which is what leaves the inward-pointing boundary flow intact.
  vars.set_rate(state_idx_log_area_sapwood_departure,
                withheld * pars.k_l - withheld_sapwood * pars.k_s +
                  rebuild_rel * (exp(-sapwood_departure) - 1.0) +
                  pars.a_sw * sapwood_marginal_);

  // [eqn 21] - Instantaneous mortality rate, now driven by relative reserves r.
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(r, vars.state(MORTALITY_INDEX)));

}

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by pars.a_y*pars.a_bio.
template <typename S>
S TF24_Strategy<S>::respiration(const S& mass_leaf, const S& mass_sapwood,
                             const S& mass_bark, const S& mass_root,
                             const S& mass_coarse_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root) +
         respiration_coarse_root(mass_coarse_root);
}

template <typename S>
S TF24_Strategy<S>::respiration_leaf(const S& mass) const {
  return pars.r_l * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_bark(const S& mass) const {
  return pars.r_b * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_sapwood(const S& mass) const {
  return pars.r_s * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_root(const S& mass) const {
  return pars.r_r * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_coarse_root(const S& mass) const {
  return pars.r_cr * mass;
}

// [eqn 14] Total turnover
template <typename S>
S TF24_Strategy<S>::turnover(const S& replacement_leaf,
                          const S& replacement_sapwood,
                          const S& mass_leaf, const S& mass_bark,
                          const S& mass_sapwood, const S& mass_root,
                          const S& mass_coarse_root) const {
   // At both fractions = 1 this is bit-identical to charging the whole flux --
   // 1.0 * x is exactly x, and the terms are summed in the same order.
   return replacement_leaf * turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          replacement_sapwood * turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root) +
          turnover_coarse_root(mass_coarse_root);
}

template <typename S>
S TF24_Strategy<S>::turnover_leaf(const S& mass) const {
  return pars.k_l * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_bark(const S& mass) const {
  return pars.k_b * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_sapwood(const S& mass) const {
  return pars.k_s * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_root(const S& mass) const {
  return pars.k_r * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_coarse_root(const S& mass) const {
  return pars.k_cr * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
template <typename S>
S TF24_Strategy<S>::net_mass_production_dt_A(const S& assimilation, const S& respiration,
                                const S& turnover) const {
  return pars.a_bio * pars.a_y * (assimilation - respiration) - turnover;
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
  // layer_thicknesses is the shared definition of dz -- do not open-code it, since
  // the vertical resistance scales with dz^2 and the two sides drifting apart
  // would be a silent squared factor neither side could detect.
  //
  // TODO(next rebaseline): pass environment.dz, the widths the environment was
  // given, instead of differencing `z`. The two agree only to rounding (on the
  // default 1.5 m over 5, layers 3-5 differ in the last bit), so the switch
  // moves every TF24 number.
  phylloptim::root_network_from_carbon(
      root_carbon_value_, phylloptim::layer_thicknesses(soil_depths_), beta_R_H,
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
  // (profit_, transpiration_, soil_consumption_, opt_psi_stem_, ...) set. Only
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
  auto optimise_at = [&](const S& PPFD_int) -> void {
    PPFD_int_used = PPFD_int;
    // Read on the energy-balance path only, and carries no row: the derivative
    // surface does not follow the energy balance.
    leaf.NIR_int_ = NIR_per_PPFD_int * odelia::util::to_passive(PPFD_int);
    // The leaf takes double, and everything the environment supplies already is
    // one. Three of its arguments are this strategy's own and carry a derivative;
    // psi_soil is the one that needs a buffer to be stripped into.
    if constexpr (std::is_same_v<S, double>) {
      leaf.set_physiology(root_network_, PPFD_int, psi_soil, soil_depths_,
                          leaf_specific_conductance_max,
                          environment.get_atm_vpd(), environment.get_ca(),
                          environment.get_leaf_temp(),
                          environment.get_atm_o2_kpa(),
                          environment.get_atm_kpa());
    } else {
      using odelia::util::to_passive;
      psi_soil_value_.resize(psi_soil.size());
      for (size_t a = 0; a < psi_soil.size(); ++a) {
        psi_soil_value_[a] = to_passive(psi_soil[a]);
      }
      leaf.set_physiology(root_network_, to_passive(PPFD_int), psi_soil_value_,
                          soil_depths_,
                          to_passive(leaf_specific_conductance_max),
                          environment.get_atm_vpd(), environment.get_ca(),
                          environment.get_leaf_temp(),
                          environment.get_atm_o2_kpa(),
                          environment.get_atm_kpa());
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
        psi_y[i]      = leaf.opt_psi_stem_ * qi;
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
      leaf.opt_psi_stem_    = function_integrator.integrate_vector(psi_y, 0.0, height);
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
  // ⚠️ GROWTH IS BILLED ON THE CARBON KEPT, `profit + shadow_cost`, NOT ON THE
  // OBJECTIVE. On TF24_floor the objective deducts `lambda_o * E`, the shadow
  // price of water: the value of the water in its best alternative use, which for
  // a leaf is assimilation later. No carbon is lost when the plant pays it; it
  // changes the aperture chosen and nothing else, which is what a Lagrange
  // multiplier does. Feeding the objective straight into growth would tax the
  // plant by carbon it never spent: measured on a 5 m plant at PPFD 1800 and theta
  // 0.25, the objective understates the carbon kept by 2.3% at lambda_o = 1e4,
  // 9.4% at 5e4, 15.4% at 1e5 and 22.6% at 2e5.
  //
  // `shadow_cost()` is exactly 0.0 on every curve but TF24_floor, and 0.0 there at
  // the default price, so this is bit-neutral at TF24's defaults. No second
  // canopy integral is needed: `shadow_cost()` reads the stored `transpiration_`,
  // which the DeepCrown branch integrated against the same weights as `profit_`,
  // and the shadow term is linear in E.
  S profit_ = leaf.profit_ + leaf.shadow_cost();
  if constexpr (!std::is_same_v<S, double>) {
    record_leaf_outputs(PPFD_int_used, psi_soil,
                        leaf_specific_conductance_max);
    profit_ = leaf_profit_ + leaf_shadow_cost_;
  }

  // --- How much would more conducting area be worth? (#516) ------------------
  // d(profit)/d(kmax) at the operating point the leaf just found.
  //
  // ⚠️ THIS COSTS TWO FURTHER LEAF SOLVES, and there is no way around it. The
  // obvious shortcut is the envelope theorem: profit is already optimised, so
  // re-evaluate at the collar potential already found and the indirect term
  // vanishes. THAT ARGUMENT IS WRONG HERE. It holds for opt_psi_stem_, which
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
    const double profit_base = leaf.profit_ + leaf.shadow_cost();
    const double dk = 1e-3 * kmax;
    psi_soil_value_.resize(psi_soil.size());
    for (size_t a = 0; a < psi_soil.size(); ++a) {
      psi_soil_value_[a] = to_passive(psi_soil[a]);
    }
    auto solve_at = [&](double conductance) {
      leaf.set_physiology(root_network_, to_passive(PPFD_int_used),
                          psi_soil_value_, soil_depths_, conductance,
                          environment.get_atm_vpd(), environment.get_ca(),
                          environment.get_leaf_temp(),
                          environment.get_atm_o2_kpa(),
                          environment.get_atm_kpa());
      solve_leaf();
    };
    solve_at(kmax + dk);
    dprofit_dkmax_ = (leaf.profit_ + leaf.shadow_cost() - profit_base) / dk;
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
      // ⚠️ THE CURVE-TYPED FORM. find_root_collar_psi() is the TF24 SHORTHAND
      // and solves TF24 whatever set_model seated, so calling it here would
      // optimise one curve, differentiate another, and report a shadow price the
      // solve never paid.
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

// [eqn 16] Fraction of production allocated to reproduction
template <typename S>
S TF24_Strategy<S>::fraction_allocation_reproduction(const S& height) const {
  return pars.a_f1 / (1.0 + exp(pars.a_f2 * (1.0 - height / pars.hmat)));
}

// Fraction of production allocated to growth
template <typename S>
S TF24_Strategy<S>::fraction_allocation_growth(const S& height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
template <typename S>
S TF24_Strategy<S>::fecundity_dt(const S& net_mass_production_dt,
                               const S& fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (pars.omega + pars.a_f3);
}

template <typename S>
S TF24_Strategy<S>::darea_leaf_dmass_live(const S& area_leaf,
                                          const S& sapwood_departure) const {
  // The coarse-root term is appended last so that at the default a_cr1 = 0 the
  // denominator is bit-identical to the pre-#349 sum.
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf, sapwood_departure)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf)
              + dmass_coarse_root_darea_leaf(area_leaf, sapwood_departure));
}

template <typename S>
S TF24_Strategy<S>::dheight_darea_leaf(const S& area_leaf) const {
  return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_leaf_darea_leaf(const S& /* area_leaf */) const {
  return pars.lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_sapwood_darea_leaf(const S& area_leaf,
                                             const S& sapwood_departure) const {
  // Priced at the Huber value the plant ACTUALLY carries. At rest
  // sapwood_per_leaf_area() returns pars.theta unchanged, so this is bitwise
  // the expression it replaces.
  return pars.rho * eta_c * pars.a_l1 * sapwood_per_leaf_area(sapwood_departure) *
    (pars.a_l2 + 1.0) * pow(area_leaf, pars.a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_bark_darea_leaf(const S& area_leaf) const {
  // At a RESTING sapwood departure deliberately: bark follows leaf area, so it
  // must not pick up the Huber value's excursion. Expressed as a_b1 times the
  // sapwood form at zero departure rather than written out, because floating
  // point multiplication is not associative -- writing the product flat changes
  // the last bits and moved a whole SCM trajectory when first tried.
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf, S(0.0));
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_root_darea_leaf(const S& /* area_leaf */) const {
  return pars.a_r1;
}

// Mass of coarse root needed for new unit area leaf, d m_cr / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_coarse_root_darea_leaf(const S& area_leaf,
                                                 const S& sapwood_departure) const {
  return pars.a_cr1 * dmass_sapwood_darea_leaf(area_leaf, sapwood_departure);
}

// Growth rate of basal diameter_stem per unit time
template <typename S>
S TF24_Strategy<S>::ddiameter_stem_darea_stem(const S& area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
template <typename S>
S TF24_Strategy<S>::area_sapwood_dt(const S& area_leaf_dt) const {
  return area_leaf_dt * pars.theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
template <typename S>
S TF24_Strategy<S>::area_heartwood_dt(const S& area_leaf,
                                      const S& sapwood_departure) const {
  return pars.k_s * area_sapwood(area_leaf, sapwood_departure);
}

// Growth rate of bark area at base per unit time
template <typename S>
S TF24_Strategy<S>::area_bark_dt(const S& area_leaf_dt) const {
  return pars.a_b1 * area_leaf_dt * pars.theta;
}

// Growth rate of stem basal area per unit time
template <typename S>
S TF24_Strategy<S>::area_stem_dt(const S& area_leaf,
                               const S& area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    // At the resting Huber value: this reporting helper is handed sizes, not a
    // state vector, so it cannot see the departure. Its callers are the R-facing
    // expansion path, which reconstructs a plant on its own trajectory.
    area_heartwood_dt(area_leaf, S(0.0));
}

// Growth rate of basal diameter_stem per unit time
template <typename S>
S TF24_Strategy<S>::diameter_stem_dt(const S& area_stem, const S& area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
template <typename S>
S TF24_Strategy<S>::mass_root_dt(const S& area_leaf,
                               const S& area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

// Growth rate of coarse-root mass per unit time
template <typename S>
S TF24_Strategy<S>::mass_coarse_root_dt(const S& area_leaf,
                                      const S& area_leaf_dt) const {
  // At the resting Huber value, for the reason area_stem_dt is.
  return area_leaf_dt * dmass_coarse_root_darea_leaf(area_leaf, S(0.0));
}

template <typename S>
S TF24_Strategy<S>::mass_live_dt(const S& fraction_allocation_reproduction,
                               const S& net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

template <typename S>
S TF24_Strategy<S>::mass_total_dt(const S& fraction_allocation_reproduction,
                                     const S& net_mass_production_dt,
                                     const S& mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

template <typename S>
S TF24_Strategy<S>::mass_above_ground_dt(const S& area_leaf,
                                       const S& fraction_allocation_reproduction,
                                       const S& net_mass_production_dt,
                                       const S& mass_heartwood_dt,
                                       const S& area_leaf_dt) const {
  const S mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  const S mass_coarse_root_dt =
    area_leaf_dt * dmass_coarse_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt - mass_coarse_root_dt;
}

template <typename S>
S TF24_Strategy<S>::mass_heartwood_dt(const S& mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


template <typename S>
S TF24_Strategy<S>::mass_live_given_height(const S& height) const {
  S area_leaf_ = area_leaf(height);
  const S mass_sapwood_ = mass_sapwood(area_sapwood(area_leaf_), height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood_ +
         mass_root(area_leaf_) +
         mass_coarse_root(mass_sapwood_);
}

template <typename S>
S TF24_Strategy<S>::height_given_mass_leaf(const S& mass_leaf) const {
  return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
}

template <typename S>
S TF24_Strategy<S>::mortality_dt(const S& relative_reserves,
                              const S& cumulative_mortality) const {

  // Growth-dependent mortality is now driven by relative NSC reserves
  // r = S/S_max (in [0,1]) rather than instantaneous productivity, so the rate
  // is bounded (see mortality_storage_dependent_dt): death becomes gradual as
  // reserves deplete instead of spiking to ~1e32 under carbon deficit (#550).
  // ⚠️ THE TEST IS THE CEILING, NOT FINITENESS, and the two must not be swapped
  // back. A cohort that failed to establish has survival exactly
  // zero, so parking its rate there is the same statement `!is_finite` used to
  // make -- and keeping the rate parked is what makes a finite hazard
  // bit-identical to the +Inf it replaces, in the state, in the rate and so in
  // the step the controller chooses. Finiteness is still tested, because a
  // hazard that arrives non-finite by any other route must park too.
  if (util::is_finite(cumulative_mortality) &&
      cumulative_mortality < establishment_failure_hazard) {
    return
      mortality_growth_independent_dt() +
      mortality_storage_dependent_dt(relative_reserves);
 } else {
    // Mortality probability is 1, so the rate calculations have nothing left to
    // describe and the state does not move again.
    return 0.0;
  }
}

template <typename S>
S TF24_Strategy<S>::mortality_growth_independent_dt() const {
  return pars.d_I;
}

// Storage-dependent growth mortality (#517), following Stefaniak et al. 2026
// (Eq 6). Bounded in [a_dG1*exp(-a_dG2), a_dG1] for r in [0,1]: full reserves
// (r=1) give near-zero excess mortality; empty reserves (r=0) give the finite
// maximum a_dG1. This boundedness is what removes the #550 ODE overflow.
template <typename S>
S TF24_Strategy<S>::mortality_storage_dependent_dt(const S& relative_reserves) const {
  return pars.a_dG1 * exp(-pars.a_dG2 * relative_reserves);
}

// NSC storage capacity: scales with sapwood mass (per Daniel, #517). mass_sapwood
// = area_sapwood(area_leaf, departure) * height * eta_c * rho.
template <typename S>
S TF24_Strategy<S>::storage_capacity(const S& area_leaf_, const S& height,
                                     const S& sapwood_departure) const {
  return pars.a_st1 *
    mass_sapwood(area_sapwood(area_leaf_, sapwood_departure), height);
}

// Initialise the storage state for a newly germinated individual at a_st3 fraction of
// its capacity (Stefaniak et al. 2026, Eq 8), so seedlings are born with
// reserves rather than starting empty (which would kill them immediately).
template <typename S>
void TF24_Strategy<S>::set_initial_states(const TF24_Environment<S>& environment,
                                       Internals<S>& vars) {
  // The seed's height is written here, not inherited from the height an
  // Individual was constructed with. Construction runs in plain arithmetic and an
  // active strategy receives its results rather than deriving them, so a height
  // that arrives that way carries no trait derivative and every rate the newborn
  // goes on to have inherits the loss. Declaring it by its own residual at this
  // point puts it where the parameters are already differentiable inputs, and the
  // value written is the same one either way.
  const SeedGeometry seed = seed_geometry();
  vars.set_state(HEIGHT_INDEX, seed.height);
  // A height written into vars leaves the slots it determines holding what the
  // constructor's plain height derived, which is the same number carrying no
  // derivative. Re-deriving them costs one allometry and is the only thing that
  // puts the seed's row into the leaf area every rate at birth size is scaled by.
  update_dependent_aux(HEIGHT_INDEX, vars);
  // A seed is on its allometry by definition, so both departures rest at zero.
  vars.set_state(state_idx_log_area_leaf_departure, S(0.0));
  vars.set_state(state_idx_log_area_sapwood_departure, S(0.0));
  vars.set_state(state_idx_storage,
                 pars.a_st3 * storage_capacity(seed.area_leaf, seed.height,
                                               S(0.0)));

  // Seed the tracked marginal balance from the plant's own first evaluation.
  // Left at zero it would start every seedling exactly at the gate's centre,
  // withholding half its replacement from birth; and a running mean has to
  // start somewhere its own history would have put it.
  //
  // Only where the gate is open. At a_pl0 = 0 the tracked value reaches no
  // rate, so the leaf solve this costs per birth would buy nothing -- and the
  // seed's carbon at birth is read from the first rate evaluation's aux anyway.
  S tracked = 0.0;
  if (odelia::util::to_passive(pars.a_pl0) > 0.0) {
    net_mass_production_dt(environment, vars);
    tracked = marginal_leaf_return_;
  }
  vars.set_state(state_idx_leaf_marginal_tracked, tracked);
}

// [eqn 20] Survival of seedlings during establishment
template <typename S>
S TF24_Strategy<S>::establishment_probability(const TF24_Environment<S>& environment) {
  const SeedGeometry seed = seed_geometry();
  return establishment_probability(
    environment,
    net_mass_production_dt(environment, seed.height, seed.area_leaf,
                           1.0 / seed.height,
                           // A germinating seed is on its allometry by
                           // definition -- height_0 is solved from seed mass,
                           // and area_leaf_0 from height_0 -- so both of its
                           // departures are the resting ones. It has no history
                           // yet, so the tracked balance starts at break-even.
                           S(0.0), S(0.0), S(0.0)));
}

// Both forms above end here. The carbon is birth-size carbon either way, whatever
// height the caller's plant happens to be at.
template <typename S>
S TF24_Strategy<S>::establishment_probability(const TF24_Environment<S>& environment,
                                               const S& net_mass_production_dt_) {

  S decay_over_time = exp(-pars.recruitment_decay * environment.time);

  if (net_mass_production_dt_ > 0) {
    const S tmp = pars.a_d0 * seed_geometry().area_leaf / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

template <typename S>
S TF24_Strategy<S>::compute_competition(const S& z, const S& height) const {
  return pars.k_I * area_leaf(height) * canopy_shape.Q_from_height(z, height);
}

// Ratio-first hot-path overload (see header): receives the cached
// competition_effect (= area_leaf(height)) and height_inverse (= 1/height), so the
// per-call area_leaf() evaluation and z/height division are hoisted out of the
// inner competition loop.
template <typename S>
S TF24_Strategy<S>::compute_competition(const S& z, const S& area_leaf_,
                                          const S& height_inverse) const {
  return pars.k_I * area_leaf_ * canopy_shape.leaf_area_above(z * height_inverse);
}

// [eqn 10] Cumulative fraction of a quantity distributed over an extent with
//          shape exponent 'eta_x', above coordinate 'z' of a total 'height'.
//          Serves the root mass distribution over soil depth.
template <typename S>
S TF24_Strategy<S>::Q(const S& z, const S& rooting_depth, const S& eta_x) const {
  if (z > rooting_depth) {
    return S(0.0);
  }
  // u^eta_x. On double the plain pow; on an active scalar the recorded eta_x
  // derivative u^eta_x * log(u) is 0 * (-inf) -- a NaN -- at u = 0, where the
  // cumulative fraction is 1, so the guard supplies that value outright.
  S u_eta;
  if constexpr (std::is_same_v<S, double>) {
    u_eta = pow(z / rooting_depth, eta_x);
  } else {
    const S u = z / rooting_depth;
    u_eta = odelia::util::to_passive(u) <= 0.0 ? S(0.0) : pow(u, eta_x);
  }
  const S tmp = 1.0 - u_eta;
  return tmp * tmp;
}

// The aim is to find a plant height that gives the correct seed mass.
template <typename S>
double TF24_Strategy<S>::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const S
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(pars.omega);

  const double tol = this->control.offspring_production_tol;
  const size_t max_iterations = this->control.offspring_production_iterations;

  auto target = [&] (double x) mutable -> S {
    return mass_live_given_height(x) - pars.omega;
  };

  if constexpr (std::is_same_v<S, double>) {
    return util::uniroot(target, h0, h1, tol, max_iterations);
  } else {
    // Bisection is affine in its bracket and blind to the residual's values, so
    // recording the search would return d(h1)/d(trait) rather than the height at
    // which mass_live equals the seed mass. Declare it by the residual instead.
    static_assert(std::is_same_v<S, double>,
                  "height_seed() finds its root by iteration; an active scalar must "
                  "declare it through implicit_value on the residual.");
  }
}

template <typename S>
void TF24_Strategy<S>::prepare_strategy() {

  // Set up the function_integrator
  function_integrator = quadrature::QK(
      // Gauss-Kronrod quadrature integeration rule (see qkrules)
      this->control.function_integration_rule);

  // Resolve the crown shading model once. The empty Control default maps to
  // TF24's own default (mean-light, its established behaviour); PPA is an
  // FF16-only stepped-light model and is rejected here.
  shading_model_ =
    shading_model_from_string(this->control.shading_model, ShadingModel::MeanLight);
  // PPA and the flat-top-box variants are FF16-only (they reshape the FF16
  // competition / light profile, which TF24 does not use).
  if (shading_model_ == ShadingModel::PPA ||
      shading_model_ == ShadingModel::FlatTopBox ||
      shading_model_ == ShadingModel::FlatTopSoftBox) {
    throw std::invalid_argument(
      "shading_model '" + this->control.shading_model +
      "' is not supported for the TF24 strategy");
  }

  // ⚠️ Sapwood acclimation needs d(profit)/d(kmax), which is measured by
  // re-solving the leaf at ONE light. Under deep-crown the profit is an integral
  // over crown positions instead, so that single-point difference is not its
  // derivative. Refuse rather than skip: skipping leaves the controller with
  // only its (negative) cost term, so the stem would shrink without bound and
  // the run would look plausible the whole way down.
  if (odelia::util::to_passive(pars.a_sw) > 0.0 &&
      shading_model_ == ShadingModel::DeepCrown) {
    throw std::invalid_argument(
      "TF24 sapwood acclimation (a_sw > 0) is not supported with the "
      "deep-crown shading model: the sapwood sensitivity is measured at a "
      "single light and deep-crown integrates over the crown");
  }

  canopy_shape.initialise(pars.eta, shading_model_);
  crown = canopy_shape.unit_crown(function_integrator);

  eta_c = CanopyShape<S>::eta_c(pars.eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
  // The residual seed_geometry closes is five composed allometries, so its slope
  // is not one term to write down: it is a tangent through the same expression,
  // taken here where the choice is visible and the scalar is already double.
  {
    using tangent = odelia::ode::tangent_scalar<double>;
    const TF24_Strategy<tangent> at_tangent = rebind<tangent>();
    tangent probe = height_0;
    odelia::ode::set_direction(probe, 1.0);
    dmass_dheight_0 = odelia::ode::derivative_along(
        at_tangent.mass_live_given_height(probe));
  }

  // theta_c is declared but NOT YET USABLE. It profiles theta along the flow
  // path, and theta is not a hydraulics-only trait: it also sets area_sapwood,
  // area_bark, their growth rates, mass_sapwood (hence construction cost,
  // respiration, turnover and NSC capacity) and the hard-coded
  // dmass_sapwood_darea_leaf derivative. Those all still read a flat pars.theta.
  //
  // Applying the profile to the hydraulic term alone would give a plant whose
  // stem conducts as though theta varied while it is built and respired as
  // though theta were constant -- two different plants sharing one trait. There
  // is no staged version of this worth having, so it is refused rather than
  // half-applied. Lift the guard in the same change that profiles theta
  // everywhere.
  //
  // Read through to_passive because prepare_strategy also runs on a rebound
  // strategy carrying an active scalar, where the comparison is of the value
  // either way and an active operand only records a statement to reach it.
  if (odelia::util::to_passive(pars.theta_c) != 0.0) {
    throw std::invalid_argument(
      "theta_c is not implemented yet: theta also sets sapwood and bark area, "
      "construction cost, respiration and storage capacity, and those still "
      "use a constant theta. A hydraulics-only theta profile would be "
      "physically inconsistent, so it is refused rather than half-applied. Use "
      "D_c to vary the height dependence of resistance.");
  }

  // A negative offset can zero the storage rate's divisor or flip its sign.
  if (!(odelia::util::to_passive(pars.storage_relaxation_offset) >= 0.0)) {
    throw std::invalid_argument("storage_relaxation_offset must be >= 0");
  }

  if (odelia::util::to_passive(stem_path_exponent()) != 0.0) {
    // L_tip is the anchor of both profiles, so it cannot be zero once either is
    // active: k_s(L) = K_s*(L/L_tip)^(2*D_c) diverges everywhere as L_tip -> 0,
    // giving zero resistance. That is a degenerate configuration to reject, not
    // a numerical edge case to tolerate.
    //
    // Written to reject zero and NaN as well as negatives: `L_tip < 0.0` alone
    // would let a zero through, which is the case this guard exists for.
    const double L_tip_value = odelia::util::to_passive(pars.L_tip);
    if (L_tip_value <= 0.0 || std::isnan(L_tip_value)) {
      throw std::invalid_argument(
        "L_tip must be > 0 when D_c or theta_c is non-zero: the within-plant "
        "profiles are defined relative to the terminal segment, and L_tip -> 0 "
        "sends sapwood-specific conductivity to infinity everywhere");
    }
    // A plant cannot be shorter than one terminal segment. Worth catching here
    // rather than downstream: height_0 is solved from seed mass, so a user
    // sweeping omega down at a fixed L_tip will eventually cross this, and the
    // symptom is a negative path length, hence a negative conductance, hence an
    // unattributable NaN twenty frames inside the leaf solver.
    if (!(L_tip_value < odelia::util::to_passive(height_0 * eta_c))) {
      throw std::invalid_argument(
        "L_tip must be shorter than the birth-size flow path (height_0*eta_c): "
        "a plant cannot be smaller than one terminal segment");
    }
  }

  if (this->is_variable_birth_rate) {
    this->extrinsic_drivers.set_variable("birth_rate", this->birth_rate_x, this->birth_rate_y,
                                         odelia::drivers::Slopes::monotone);
  } else {
    this->extrinsic_drivers.set_constant("birth_rate", this->birth_rate_y[0]);
    birth_rate_ = this->birth_rate_y[0];
  }
  if constexpr (std::is_same_v<S, double>) {
    // (P50, c) per curve, and phylloptim derives b and psi_crit itself. Handing
    // over the derived pair instead would state the curve twice and let the two
    // disagree silently, which is why they are not in the trait set.
    leaf = Leaf(pars.vcmax_25, pars.stem_c, pars.stem_P50,
                pars.root_c, pars.root_P50,
                pars.TF24_beta2, pars.jmax_25, pars.phi_J,
                pars.curv_fact_elec_trans, pars.curv_fact_colim,
                this->control.GSS_tol_abs, this->control.vulnerability_curve_ncontrol,
                this->control.ci_abs_tol, this->control.ci_niter,
                pars.TF24_cost_scale);
    // Not a constructor argument, and written before any physiology is set, so
    // the temperature block derives R_d_ from it on the first call.
    leaf.R_d_25 = pars.R_d_25;
    // Penman-Monteith leaf energy balance (#523): enable per pars (default off,
    // backward-compatible) and pass the leaf-dimension trait. Wind speed is a
    // per-timestep driver, set from the environment before each set_physiology.
    //
    // ⚠️ BOTH ARE SETTABLE AND REACH NOTHING IF THIS IS DROPPED. Moving the model
    // body out of src/tf24_strategy.cpp and into this header lost the pair, and
    // the symptom was Tleaf reporting the air temperature exactly on every run
    // with the balance switched on -- an answer, not an error.
    leaf.use_energy_balance_ = (odelia::util::to_passive(pars.use_energy_balance) != 0.0);
    leaf.d_ = odelia::util::to_passive(pars.d);
    leaf.alpha_PAR_ = odelia::util::to_passive(pars.alpha_PAR);
    leaf.alpha_NIR_ = odelia::util::to_passive(pars.alpha_NIR);
    // ⚠️ SEATED ON TF24_floor, WHICH IS TF24 AT lambda_o = 0 -- and exactly, not
    // approximately: each term is TF24's own expression, so zeroing the price
    // adds an exact zero to TF24's exact value. The default price IS zero, so no
    // run moves; what changes is that the price now reaches the model at all.
    // Without this the parameter is settable and inert, which is the silent
    // failure the whole (P50, c) reparameterisation was about avoiding.
    leaf.TF24_floor_lambda_o = pars.TF24_floor_lambda_o;
    leaf.set_model(phylloptim::Leaf::CostCurve::TF24_floor, true);
  } else {
    static_assert(std::is_same_v<S, double>,
                  "Leaf carries double; an active strategy must supply the leaf's "
                  "local Jacobian across this boundary, not template Leaf.");
  }
}

template <typename S>
typename TF24_Strategy<S>::ptr make_strategy_ptr(TF24_Strategy<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24_Strategy<S> >(s);
}

}

#endif
