// Built from  inst/include/plant/models/ff16_strategy.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_STRATEGY_H_
#define PLANT_PLANT_TF24_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/qag.h>
#include <plant/leaf_model.h>
#include <plant/canopy_shape.h>
#include <plant/stem_hydraulics.h>

namespace plant {

// Biological (user-settable) parameters for the TF24 strategy. Held as a value
// member `pars` on TF24_Strategy and exposed to R as a nested RcppR6 list class
// (access as `s$pars$lma`). Only parameters that were previously exposed to R
// live here; derived quantities (eta_c, height_0, ...), the embedded Leaf
// model, solver tolerances and hard-coded hydraulic-root constants stay as
// plain members on the strategy.
struct TF24_Pars {
  // * Core traits
  double lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  double rho       = 608.0;      // Wood density [kg/m3]
  double hmat      = 16.5958691; // Height at maturation [m]
  double omega     = 3.8e-5;     // Seed mass [kg]
  // * Individual allometry
  double eta       = 12.0;       // Canopy shape parameter [dimensionless]
  double theta     = 1.0/4669;   // Sapwood area per leaf area [dimensionless]
  double a_l1      = 5.44;       // height with 1m2 leaf [m]
  double a_l2      = 0.306;      // scaling of height with leaf area
  double a_r1      = 0.07;       // Root mass per leaf area [kg / m]
  double a_b1      = 0.17;       // Ratio of bark area : sapwood area
  // * Production
  double r_s    = 4012.0 / 608.0; // Sapwood respiration per stem mass
  double r_b    = 2.0 * r_s;      // Bark respiration (assumed 2 x sapwood)
  double r_r    = 217.0;          // Root respiration per mass
  double r_l    = 39.27 / 0.1978791; // Leaf dark respiration per leaf mass
  double a_y    = 0.7;            // Carbon conversion parameter
  double a_bio  = 2.45e-2;        // CO2 -> dry mass [kg / mol]
  double k_l    = 0.4565855;      // Leaf turnover [/yr]
  double k_b    = 0.2;            // Bark turnover [/yr]
  double k_s    = 0.2;            // Sapwood turnover [/yr]
  double k_r    = 1.0;            // Root turnover [/yr]
  double a_p1   = 151.177775377968;   // LRC hyperbola [mol CO2 / yr / m2]
  double a_p2   = 0.204716166503633;  // LRC hyperbola shape
  // * Seed production
  double a_f3   = 3.0 *  3.8e-5;  // Accessory cost of reproduction [kg/seed]
  double a_f1   = 1.0;            // Maximum allocation to reproduction
  double a_f2   = 50;             // Size range across which individuals mature
  // * Mortality parameters
  double S_D    = 0.25;           // Probability of survival during dispersal
  double a_d0   = 0.1;            // Parameter for seedling survival
  double d_I    = 0.01;           // Baseline intrinsic mortality [/yr]
  // a_dG1 / a_dG2 now shape the *storage-dependent* growth mortality (#517):
  //   mortality_storage_dependent_dt(r) = a_dG1 * exp(-a_dG2 * r),  r = S/S_max.
  // Previously these fed the instantaneous productivity mortality
  // a_dG1*exp(-a_dG2*productivity_area), which overflowed to ~1e32 under deep
  // carbon deficit (the #550 blow-up). Buffering through the NSC pool bounds the
  // argument to r in [0,1], so this term is now bounded in [a_dG1*e^-a_dG2, a_dG1].
  double a_dG1  = 5.5;            // Max growth-related (low-reserve) mortality [/yr]
  double a_dG2  = 20.0;           // Sensitivity of mortality to relative reserves
  // * NSC storage pool (#517) -- buffers growth & mortality against short-term
  //   productivity swings. See Stefaniak et al. 2026 (plantNSC) for the design.
  //   Carbon first charges storage; growth/reproduction are mobilised from it,
  //   but gated so growth only proceeds when reserves are ample (a plant should
  //   not grow down its stores). Reserves-based mortality rises as they deplete.
  double a_st1  = 0.10;           // Storage capacity per unit sapwood mass [kg NSC / kg]
  double a_st2  = 0.1;            // Reserve fraction at which growth is half-on [0-1]
  double a_st3  = 0.8;            // Initial storage at birth [fraction of capacity]
  // * Light capture
  double k_I = 0.5;
  // * Leaf hydraulic / photosynthesis traits (default Eucalyptus saligna)
  double vcmax_25 = 96;
  double p_50 = 1.85;
  // Sapwood-specific conductivity of the TERMINAL segment. Was 1, a whole-stem
  // value under the height-linear model; back-derived to the tip by
  // TF24_K_s_from_whole_stem(1) so that resistance is UNCHANGED at
  // TF24_H_ANCHOR = 1 m. The ratio is 2.9782.
  //
  // Only theta/K_s is identifiable, so agreement with the height-linear model
  // holds at exactly ONE height and nowhere else. That is by construction:
  // changing the height dependence is the object of the exercise.
  double K_s = 0.33577377016801868;
  double c = log(log(1-0.5)/log(1-0.88))/(log(p_50) - log(5.16));
  double b = p_50 / std::pow(-log(1 - 50.0 / 100.0), 1 / c);
  double psi_crit = b*std::pow(log(1/0.05),1/c); // derived from b and c
  double beta1 = 20000;
  double beta2 = 1.5;
  double g1_TF24 = 7.5;
  double jmax_25 = vcmax_25*1.64;
  double a = 0.30; // effective quantum yield of electron transport
  double curv_fact_elec_trans = 0.7;
  double curv_fact_colim = 0.99;
  double var_sapwood_volume_cost = 1;
  // nitrogen allocation traits (parameterised from Austraits 4.1.0)
  double nmass_l = 13e-3; // kg N kg^-1 mass
  double nmass_s = 1.98e-3; // kg N kg^-1 mass
  double nmass_b = 3.40e-3; // kg N kg^-1 mass
  double nmass_r = 3.35e-3; // kg N kg^-1 mass
  double dmass_dN = 0; // change in mass per change in kg kg^-1 N
  // shape exponent for the Q() root-fraction-with-depth profile
  double root_depth_shape_eta = 0.2;
  // * Root hydraulics
  // Root vulnerability curve, proportion of conductivity =
  // exp(-(psi/root_b)^root_c). Declared before root_psi_crit, which derives
  // from them. Previously fixed members of TF24_Strategy and unreachable from
  // R, which pinned root shutoff at ~5.87 MPa -- too conservative for taxa that
  // operate below that (e.g. Acacia aneura), and unavailable for calibration.
  double root_c = 2.680147;
  double root_b = 3.898245;
  // Potential at 5% remaining root conductivity [MPa]. Derived, exactly as
  // psi_crit is from b and c: if you set root_b or root_c directly, set this
  // too, or the vulnerability curve and the shutoff threshold disagree.
  double root_psi_crit = root_b*std::pow(log(1/0.05),1/root_c);
  // Maximum rooting depth [m]. Rooting depth is min(height, rooting_depth_max),
  // so this also bounds the depth over which roots can draw water. Should not
  // exceed the soil column depth (TF24_Environment `depth`, default 1.5 m):
  // layers below the column do not exist, so deepening roots alone gains
  // nothing without deepening the soil as well.
  double rooting_depth_max = 1.5;
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
  double D_c = 0.2;
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
  double theta_c = 0.0;
  // Terminal segment length [m]. Not an innocuous numerical cutoff: it enters
  // the resistance with elasticity beta ~ 0.6 and trades off exactly against
  // both K_s and theta, which are identifiable only in the grouping
  // theta*L_tip^beta/K_s. It must therefore come from the SAME terminal-segment
  // definition over which K_s and theta were measured -- none of the three may
  // be calibrated independently of the others (invariance criterion I5).
  double L_tip = 0.02;
  // Germination
  double recruitment_decay = 0.0;
  // Penman-Monteith leaf energy balance (#523). use_energy_balance gates PM
  // (0 = off, today's Tleaf=Tair behaviour; != 0 = on); default off preserves
  // backward compatibility. d is the characteristic leaf dimension (m) for the
  // aerodynamic resistance ra = C_ra*sqrt(d/U0); inert while PM is off.
  double use_energy_balance = 0.0;
  double d = 0.05;
};

class TF24_Strategy: public Strategy<TF24_Environment> {
public:
  typedef std::shared_ptr<TF24_Strategy> ptr;
  TF24_Strategy();

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
  // v9 (#615): stem resistance becomes a path integral over two within-plant
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
  // which reproduces this file's pre-v9 pinned values unmodified.
  //
  // THE SIGN OF THE EFFECT DEPENDS ON DENSITY, which is the most important
  // thing to know about this change. Individually, plants are better off
  // wherever they are taller than the anchor (assimilation ratio 0.997 at 0.5 m,
  // 1.001 at 1 m, 1.055 at 5 m, 1.143 at 10 m, single plant, wet soil, no
  // competition). But lower resistance also means faster transpiration, so in a
  // dense stand everyone draws the shared soil column down faster and the patch
  // does WORSE. One-species SCM, hmat = 5, max_patch_lifetime = 5:
  //
  //     birth_rate    0.5      2       20
  //     ratio        1.216   1.095   0.857
  //
  // The pinned scenarios below all run at birth_rate = 20, i.e. at the least
  // favourable end of that range; they are not representative of the change's
  // sign in general.
  //
  //     one-species SCM offspring     81.8571 -> 70.1963    -14.25%
  //     two-species, fast             67.3220 -> 56.2809    -16.40%
  //     two-species, slow           2.7745e-4 -> 1.2629e-4  -54.49%
  //     birth-date coordinate, fast  287.1604 -> 273.8933    -4.62%
  //     birth-date coordinate, slow   59.5320 ->  49.3474   -17.11%
  //     seeded stochastic counts        81 / 3 -> 104 / 3
  //
  // The hydraulic gateway runs longer patches at the default hmat and moves the
  // other way, up by 2.2x to 320x on every scenario, with S02 crossing R0 = 1 so
  // persistence goes 2/8 -> 3/8. 8/8 still run, 0 crash, no `observed` flips.
  //
  // theta_c is declared but REFUSED (prepare_strategy throws on any non-zero
  // value). theta is read by the carbon budget as well as the hydraulic term,
  // so a hydraulics-only profile would be an incoherent model rather than a
  // staging step; it lands everywhere at once or not at all.
  //
  // p_50 is deliberately UNCHANGED (2.8887 MPa): make_TF24_hyperpar derives the
  // vulnerability curve from K_s, so B_Hv1 was re-anchored 0.4607063 ->
  // 0.36591565341924093 to stop the reparameterisation from also moving it.
  // Left alone it would have gone to 3.5933 MPa.
  static constexpr int scientific_version = 9;

  double compute_average_light_environment(double z, double height,
                                           const TF24_Environment &environment);

  // calculate the amount of water transpired relativised by leaf area index.

  double evapotranspiration_dt(double area_leaf_, int soil_layer);


  // Overrides ----------------------------------------------

  // update this when the length of state_names changes
  static size_t state_size () { return 6; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
      "area_heartwood",
      "mass_heartwood",
      "storage"
      });
  }

  std::vector<std::string> aux_names() {
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
      "stom_cond_CO2",
      // Net CO2 assimilation at the optimal operating point, per unit leaf
      // area (umol CO2 m^-2 s^-1). Net, not gross: Leaf::assim_colimited()
      // subtracts dark respiration R_d_, so gross = assimilation + R_d_ with
      // R_d_ = 0.015 * vcmax_ at the acclimated vcmax_.
      "assimilation"
    });
    // add the associated computation to compute_rates and compute there
    if (collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  // Translate generic methods to TF24 strategy leaf area methods

  double competition_effect(double height) const {
    return area_leaf(height);
  }

  void refresh_indices();


  // TF24 Methods  ----------------------------------------------

  // [eqn 2] area_leaf (inverse of [eqn 3])
  double area_leaf(double height) const;

  // [eqn 1] mass_leaf (inverse of [eqn 2])
  double mass_leaf(double area_leaf) const;

  // [eqn 4] area and mass of sapwood
  double area_sapwood(double area_leaf) const;
  double mass_sapwood(double area_sapwood, double height) const;

  // [eqn 5] area and mass of bark
  double area_bark(double area_leaf) const;
  double mass_bark (double area_bark, double height) const;

  double area_stem(double area_bark, double area_sapwood,
                            double area_heartwood) const;
  double diameter_stem(double area_stem) const;

  // [eqn 7] Mass of (fine) roots
  double mass_root(double area_leaf) const;

  // [eqn 8] Total Mass
  double mass_live(double mass_leaf, double mass_bark,
                   double mass_sapwood, double mass_root) const;

  double mass_total(double mass_leaf, double mass_bark, double mass_sapwood,
                    double mass_heartwood, double mass_root) const;

  // Above-ground mass = leaf + all stem components (bark + sapwood +
  // heartwood); excludes roots.
  double mass_above_ground(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_heartwood) const;

  void compute_rates(const TF24_Environment& environment,
                Internals& vars);
  
  void compute_roots(const TF24_Environment& environment,
                Internals& vars);

  void update_dependent_aux(const int index, Internals& vars);

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  double assimilation(const TF24_Environment& environment, double height,
                      double area_leaf);
  // [Appendix S6] Per-leaf photosynthetic rate.
  double assimilation_leaf(double x) const;

  // [eqn 13] Total maintenance respiration
  double respiration(double mass_leaf, double mass_sapwood,
                     double mass_bark, double mass_root) const;

  double respiration_leaf(double mass) const;
  double respiration_bark(double mass) const;
  double respiration_sapwood(double mass) const;
  double respiration_root(double mass) const;

  // [eqn 14] Total turnover
  double turnover(double mass_leaf, double mass_bark,
                  double mass_sapwood, double mass_root) const;
  double turnover_leaf(double mass) const;
  double turnover_bark(double mass) const;
  double turnover_sapwood(double mass) const;
  double turnover_root(double mass) const;

  // [eqn 15] Net production
  double net_mass_production_dt_A(double assimilation, double respiration,
                                  double turnover) const;

  virtual double net_mass_production_dt(const TF24_Environment& environment,
                                double height, double area_leaf_,
                                double height_inverse);

  // Resolve the leaf operating point on the already-set-up `leaf` (i.e. after
  // leaf.set_physiology(...)). Base TF24 optimises the root-collar psi via
  // golden-section search; the TF24f variant overrides this to make the optimum
  // chase a tracked ODE state (#525). Called per crown light point from
  // net_mass_production_dt, so it must be virtual to dispatch to the override
  // when net_mass_production_dt is reused unchanged by the subclass.
  virtual void solve_leaf();
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // height state and the cached aux slots itself, so the generic Individual
  // does not need to know TF24's state/aux layout.
  double net_mass_production_dt(const TF24_Environment& environment,
                                const Internals& vars) {
    return net_mass_production_dt(environment, vars.state(HEIGHT_INDEX),
                                  vars.aux(aux_idx_competition_effect),
                                  vars.aux(aux_idx_height_inverse));
  }

  // [eqn 16] Fraction of whole plan growth that is leaf
  virtual double fraction_allocation_reproduction(double height) const;
  double fraction_allocation_growth(double height) const;
  // [eqn 17] Rate of offspring production
  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  double darea_leaf_dmass_live(double area_leaf) const;

  // change in height per change in leaf area
  double dheight_darea_leaf(double area_leaf) const;
  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  double dmass_leaf_darea_leaf(double area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  double dmass_sapwood_darea_leaf(double area_leaf) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  double dmass_bark_darea_leaf(double area_leaf) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  double dmass_root_darea_leaf(double area_leaf) const;
  // Growth rate of basal diameter_stem per unit stem area
  double ddiameter_stem_darea_stem(double area_stem) const;
  // Growth rate of components per unit time:
  double area_leaf_dt(double area_leaf_dt) const;
  double area_sapwood_dt(double area_leaf_dt) const;
  double area_heartwood_dt(double area_leaf) const;
  double area_bark_dt(double area_leaf_dt) const;
  double area_stem_dt(double area_leaf, double area_leaf_dt) const;
  double diameter_stem_dt(double area_stem, double area_stem_dt) const;
  double mass_root_dt(double area_leaf,
                       double area_leaf_dt) const;
  double mass_live_dt(double fraction_allocation_reproduction,
                       double net_mass_production_dt) const;
  double mass_total_dt(double fraction_allocation_reproduction,
                        double net_mass_production_dt,
                        double mass_heartwood_dt) const;
  double mass_above_ground_dt(double area_leaf,
                               double fraction_allocation_reproduction,
                               double net_mass_production_dt,
                               double mass_heartwood_dt,
                               double area_leaf_dt) const;

  double mass_heartwood_dt(double mass_sapwood) const;

  double mass_live_given_height(double height) const;
  double height_given_mass_leaf(double mass_leaf_) const;


  double mortality_dt(double relative_reserves, double cumulative_mortality) const;
  double mortality_growth_independent_dt()const ;
  // Storage-dependent growth mortality (#517): rises smoothly as relative
  // reserves r = S/S_max deplete, bounded in [a_dG1*e^-a_dG2, a_dG1].
  double mortality_storage_dependent_dt(double relative_reserves) const;
  // NSC storage capacity S_max = a_st1 * mass_sapwood [kg NSC].
  double storage_capacity(double area_leaf, double height) const;
  // Seed the storage state for a newly germinated individual (#517).
  void set_initial_states(const TF24_Environment& environment, Internals& vars);
  // [eqn 20] Survival of seedlings during establishment, from the carbon a
  // seedling produces at birth size. This form works that carbon out.
  double establishment_probability(const TF24_Environment& environment);
  // The same, for a newborn whose rates have just been computed. A newborn is
  // already at birth size, so compute_rates has left that carbon in aux and the
  // leaf need not be solved there twice.
  double establishment_probability(const TF24_Environment& environment,
                                   const Internals& vars) {
    return establishment_probability(environment,
                                     vars.aux(aux_idx_net_mass_production_dt));
  }
  // The equation the two above share.
  double establishment_probability(const TF24_Environment& environment,
                                   double net_mass_production_dt_);

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  double compute_competition(double z, double height) const;
  // Optimised overload called from Individual<TF24>::compute_competition with the
  // cached competition_effect (= area_leaf(height)) and height_inverse (= 1/height)
  // aux values, matching the shared individual.h interface (no recompute per call).
  double compute_competition(double z, double area_leaf_,
                             double height_inverse) const;
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  double compute_competition(double z, const Internals& vars) const {
    return compute_competition(z, vars.aux(aux_idx_competition_effect),
                               vars.aux(aux_idx_height_inverse));
  }

  // The fraction of root mass below soil depth `z`, for a plant rooted to
  // `rooting_depth` with shape exponent `eta_x` (pars.root_depth_shape_eta). The
  // canopy's own cumulative form is CanopyShape::Q, at pars.eta.
  double Q(double z, double rooting_depth, double eta_x) const;

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // Set constants within TF24_Strategy
  void prepare_strategy();

  // Birth height of a (germinated) seed. Strategy-agnostic accessor used by
  // the templated Individual; here height_0 is derived in prepare_strategy().
  double initial_height() const { return height_0; }

  // Crown shading model, resolved once from control.shading_model in
  // prepare_strategy(). TF24 supports deep-crown, mean-light (its default)
  // and crown-centre; PPA is not available for TF24.
  ShadingModel shading_model_ = ShadingModel::MeanLight;

  // Biological (user-settable) parameters; see TF24_Pars above.
  TF24_Pars pars;

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  double eta_c     = NA_REAL; // crown shape factor, precomputed from pars.eta
  // beta = 2*pars.D_c + pars.theta_c, the exponent of the stem path integral:
  // leaf-specific resistance grows as H^(1-beta). Named for the path rather than
  // for either parameter, since it belongs to neither. Cached because
  // prepare_strategy() is also where the (beta, L_tip) precondition is checked,
  // and the two belong together: nothing may reach
  // stem_hydraulics::effective_path_length that has not been through that
  // validation.
  double stem_path_exponent_ = 0.0;
  CanopyShape canopy_shape;
  // Height and leaf area of a (germinated) seed
  double height_0  = NA_REAL;
  double area_leaf_0;

  // Embedded leaf hydraulic/photosynthesis sub-model, built in prepare_strategy()
  Leaf leaf;

  // Width of the smooth reserve gate G(r) on growth (#517); small relative to
  // [0,1] so the switch about the growth threshold a_st2 is fairly sharp but
  // differentiable. storage_prod_eps smooths the positive-part of net production
  // (replacing the old hard net>0 growth cutoff) for AD-readiness.
  double storage_gate_width = 0.1;
  double storage_prod_eps   = 1e-4;

  // Solver tolerances and other constants not currently exposed to R
  double newton_tol_abs = 0.001;
  double GSS_tol_abs = 1e-3;
  double vulnerability_curve_ncontrol = 100;
  double ci_abs_tol = 1e-6;
  double ci_niter = 1000;
  // The root-architecture model's two constants. They used to be handed to the
  // Leaf constructor; since phylloptim #33 the leaf takes the RESISTANCES and this
  // strategy owns the model that produces them, which is where they belong -- the
  // 1/3 : 2/3 root split and the dz^2 vertical scaling were never gas exchange.
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
  int aux_idx_stom_cond_CO2 = -1;
  int aux_idx_assimilation = -1;
  int aux_idx_area_sapwood = -1;       // only present when collect_all_auxiliary
  int state_idx_area_heartwood = -1;
  int state_idx_mass_heartwood = -1;
  int state_idx_storage        = -1;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;

  // Reusable per-layer root-mass buffer, refilled (not reallocated) each
  // net_mass_production_dt call to avoid a heap allocation per derivs eval.
  std::vector<double> root_carbon_per_leaf_area_;

  // And the resistances derived from it, held the same way and for the same
  // reason (phylloptim #33). Refilled through root_network_from_carbon's in-place
  // overload: building five fresh vectors per solve measured +0.074 us there,
  // about +2% of a whole solve. It must NOT be moved from -- phylloptim's
  // set_root_network takes it by const reference precisely so this buffer keeps
  // its capacity across calls.
  phylloptim::RootNetwork root_network_;
};

TF24_Strategy::ptr make_strategy_ptr(TF24_Strategy s);

}

#endif
