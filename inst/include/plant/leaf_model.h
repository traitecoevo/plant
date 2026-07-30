// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

// TODO: replace with constants
// #define umol_per_mol_to_mol_per_mol 0.000...
// #define umol_per_mol_to_Pa ...
// #define kg_to_mol_h2o ...
// #define kPa_to_Pa ...

#include <plant/models/tf24_environment.h>
#include <plant/qag.h>
#include <plant/uniroot.h>
#include <plant/optimize.h>

namespace plant {

// double check best namespace for constants (private vs global)
// converts vcmax to jmax 25 unitless (Sperry et el. (2017))
static const double vcmax_25_to_jmax_25 = 1.67;

// kJ mol ^-1
static const double vcmax_ha = 60000;
// kJ mol ^-1
static const double vcmax_H_d = 200000;
// kJ mol ^-1
static const double vcmax_d_S = 650;

// kJ mol ^-1
static const double jmax_ha = 30000;
// kJ mol ^-1
static const double jmax_H_d = 200000;
// kJ mol ^-1
static const double jmax_d_S = 650;

// umol ^ -1 mol ^ 1
static const double gamma_25 = 42.75;
// dimensionless
static const double gamma_c = 19.02;
// kJ mol ^-1
static const double gamma_ha = 37.83e3;

// umol mol ^-1
static const double kc_25 = 404.9 ;
// dimensionless
static const double kc_c = 38.05;
// kJ mol ^-1
static const double kc_ha = 79.43e3;

// umol mol ^-1
static const double ko_25 = 278400 ;
// dimensionless
static const double ko_c = 20.30;
// kJ mol ^-1
static const double ko_ha = 36.38e3;

// Pa umol ^ -1 mol ^ 1 
static const double umol_per_mol_to_Pa = 0.1013;

// mol H2o kg ^-1
static const double kg_to_mol_h2o = 55.4939;
// kg mol^-1: molar mass of water, for converting molar water flux back to kg.
// (Intentionally distinct from 1/kg_to_mol_h2o, which it does not exactly equal;
// kept at the historical 0.018015 to preserve results.)
static const double kg_per_mol_h2o = 0.018015;
// mol mol ^-1 / (umol mol ^-1)
static const double umol_to_mol = 1e-6;
// Pa kPa^-1
static const double kPa_to_Pa = 1000.0;

// universal gas constant J mol^-1 K^-1
static const double R = 8.314;

//convert deg C to deg K
static const double C_to_K = 273.15;

//H20:CO2 stomatal diffusion ratio
static const double H2O_CO2_stom_diff_ratio = 1.67;

const double gravity_head = 9.8e-3; // MPa / m

// number of intergration steps
const double n = 5;

// --- Penman-Monteith leaf energy balance (minimal core; #523) -----------------
// See notes/penman-monteith/. These back Leaf::leaf_temp_from_E and the es/Delta
// helpers. Only used on the (default-off) use_energy_balance_ path.
// latent heat of vaporisation of water, J kg^-1 (fixed at 25 deg C)
static const double latent_heat_vap = 2.45e6;
// volumetric heat capacity of air, J m^-3 K^-1
static const double vol_heat_cap_air = 1200.0;
// PAR energy conversion: ~4.57 umol photons per J of PAR (shortwave)
static const double umol_par_per_joule = 4.57;
// shortwave absorbed ~= 2 * absorbed PAR (PAR ~= 50% of shortwave; doc 3.3)
static const double sw_abs_per_par = 2.0;
// fixed net longwave (cooling) offset, W m^-2 (clear-sky approximation; doc 3.3).
// A faithful treatment makes outgoing longwave depend on Tleaf (would make Rn
// implicit in Tleaf); staged for a sensitivity gate and tracked in #581.
static const double longwave_net_offset = -40.0;
// fixed aerodynamic resistance fallback, s m^-1 (doc 6/7.4; used when the wind
// model is unavailable, e.g. a bare Leaf with no wind/d set)
static const double aerodynamic_resistance_fixed = 50.0;
// leaf boundary-layer coefficient C_ra, s^0.5 m^-1, in ra = C_ra*sqrt(d/U) (doc 4.1)
static const double aerodynamic_resistance_coef = 200.0;
// Physical clamp on the energy-balance leaf temperature (deg C). The linear
// balance can return absurd temperatures at non-equilibrium operating points an
// optimiser may probe (a large transpiration driving Tleaf below absolute zero,
// making the Arrhenius block non-finite). Real leaves operate far inside this
// range; clamping keeps A(Tleaf) finite so such points get a finite (poor)
// profit and are simply rejected.
static const double leaf_temp_min = -40.0;
static const double leaf_temp_max = 70.0;

class Leaf {
public:
  //anonymous Leaf function as in canopy.h
  Leaf();
  
  Leaf(double vcmax_25, 
       double c, 
       double b, 
       double psi_crit,
       double root_c,
       double root_b,
       double root_psi_crit,
       double beta2, 
       double jmax_25, 
       double a, 
       double curv_fact_elec_trans, 
       double curv_fact_colim,
       double GSS_tol_abs,
       double vulnerability_curve_ncontrol,
       double ci_abs_tol,
       double ci_niter,
      double g1_TF24,
    double beta_R_H,
    double beta_R_V); 
        
  quadrature::QAG integrator;
  odelia::interpolator::Interpolator transpiration_from_psi;
  odelia::interpolator::Interpolator psi_from_transpiration;
  // pre-computed root vulnerability curve (same role as transpiration_from_psi for xylem)
  odelia::interpolator::Interpolator root_vuln_from_psi;
  // cumulative integral of the root vulnerability curve, G(m) = int_0^m f_r(s) ds,
  // indexed by magnitude m = -psi. Lets E_from_Soil_to_Root_Collar obtain the
  // mean conductivity over a potential interval from 2 evals instead of (n+1)
  // (same pre-integrated-curve trick as transpiration_from_psi for the xylem).
  odelia::interpolator::Interpolator root_vuln_integral_from_psi;

  // psi_from_E

  double vcmax_25;
  double c;
  double b;
  double psi_crit;  // derived from b and c
  double root_c;
  double root_b;
  double root_psi_crit;
  double beta2;
  double jmax_25;
  double a;
  double curv_fact_elec_trans; // unitless - obtained from Smith and Keenan (2020)
  double curv_fact_colim;
  double GSS_tol_abs;
  double vulnerability_curve_ncontrol;
  double ci_abs_tol;
  double ci_niter;
  double g1_TF24;
  double beta_R_H;
  double beta_R_V;
  double soil_number_of_depths_;
  int max_soil_layer;

  double ci_;
  double stom_cond_CO2_;
  double assim_colimited_;
  double transpiration_;
  std::vector<double> soil_consumption_;
  double profit_;
  double psi_stem;
  double lambda_;
  double lambda_analytical_;
  double hydraulic_cost_;
  
  double electron_transport_;
  double gamma_;
  double ko_;
  double kc_;
  double km_;
  double R_d_;
  double leaf_specific_conductance_max_;
  double sapwood_volume_per_leaf_area_;
  double k_s_;
  double root_mass_;
  std::vector<double> c_r_V_;
  std::vector<double> c_r_H_;
  double area_leaf_;
  double rho_;
  double vcmax_;
  double jmax_;
  double lma_; //kg m^-2
  double a_bio_;
  
  std::vector<double> psi_soil_;
  std::vector<double> psi_soil_inverted_;
  // Per-layer cache of root_vuln_integral_from_psi.eval(-psi_soil_inverted_[i]).
  // psi_soil_inverted_ is fixed for the whole find_root_collar_psi solve, so the
  // soil-side endpoint of the cumulative-integral lookup in
  // E_from_Soil_to_Root_Collar is constant across every (re)evaluation of the
  // nested root-finders. Precomputing it once per solve (alongside the
  // P_x_r-side eval, hoisted out of the layer loop) collapses ~2 spline evals
  // per layer to ~1 per call. Rebuilt in find_root_collar_psi.
  std::vector<double> root_vuln_integral_soil_;
  std::vector<double> soil_depth_;
  std::vector<double> z_soil_mid_;
  // Per-layer gravitational head gravity_head * z_soil_mid_[i], precomputed once
  // per solve in set_physiology (z_soil_mid_ is fixed across find_root_collar_psi).
  // Used three times per layer in E_from_Soil_to_Root_Collar's hot loop; caching
  // it removes a redundant multiply per layer per (re)evaluation.
  std::vector<double> grav_head_z_;
  bool use_precomputed_z_soil_mid_;
  double dz_;
  std::vector<double> r_R_H_min;
        // vertical root resistance
    std::vector<double> r_R_V;

            // cumulative vertical sum of root resistance
    std::vector<double> r_R_V_sum;
    

  double leaf_temp_;
  // Penman-Monteith leaf energy balance state (#523), only meaningful on the
  // use_energy_balance_ path. Set once per set_physiology; Tleaf itself is a
  // per-operating-point quantity computed in set_leaf_states_rates_from_psi_stem.
  double Tair_;  // air temperature, deg C (reinterprets the leaf_temp driver)
  double Rn_;    // net radiation at the leaf, W m^-2
  double ra_;    // aerodynamic (boundary-layer) resistance, s m^-1
  // Gate for the PM leaf energy balance. Default OFF: today's path runs
  // (prescribed leaf_temp, single-shot cached Arrhenius). R-settable (#523 full
  // cut) so TF24 (via pars.use_energy_balance) and the leaf-level demo can
  // turn PM on; default preserves backward compatibility.
  bool use_energy_balance_ = false;
  // Boundary-layer inputs for ra = C_ra*sqrt(d/U0) (doc 4.1). d is a per-strategy
  // trait (set from pars.d in prepare_strategy); wind_speed_ is the per-timestep
  // above-canopy driver (set from the environment before set_physiology). Both
  // R-settable so a bare Leaf can exercise the wind model; if either is unusable
  // set_physiology falls back to the fixed ra. Only read on the PM path.
  double d_ = 0.05;          // characteristic leaf dimension, m
  double wind_speed_ = 2.0;  // above-canopy wind speed U0, m s^-1
  double PPFD_;
  double atm_vpd_;
  double atm_o2_kpa_;
  double atm_kpa_;
  double ca_;
  double root_collar_psi_;
  double assim_max_;

  
  double opt_psi_stem_;
  double opt_ci_;
  double count;
  double E_up_;

  // --- Medlyn stomatal-conductance model (from develop #450) ------------------
  // Standalone, R-callable alternative to the root-collar profit optimisation
  // (solve_medlyn_ci_*); NOT used by the TF24 compute path, which optimises
  // psi_stem directly. g0/g1 default to the published values and are exposed as
  // settable fields. beta_ (soil-moisture stress) uses theta_/theta_w_/theta_fc_,
  // set from the default soil-moisture members in set_physiology.
  double g0 = 0.022;        // residual stomatal conductance (umol m^-2 s^-1)
  double g1 = 2.57;         // sensitivity to vpd (kPa^0.5)
  double medlyn_model_gs_;  // mol CO2 m^-2 s^-1
  double theta_w_;          // current soil water content at wilting point (m^3 m^-3)
  double theta_fc_;         // current soil water content at field capacity (m^3 m^-3)
  double theta_;            // current soil water content (m^3 m^-3)

  // 1-entry memo for transpiration(). Within a single root-collar/profit solve
  // leaf_specific_conductance_max_ and the transpiration spline are fixed, so
  // supply-side transpiration depends only on (psi_stem, psi_upstream). That
  // pair is queried repeatedly with identical values per profit evaluation
  // (psi_stem_to_ci -> stom_cond_CO2 -> transpiration, then transpiration again
  // for transpiration_). Caching the last result avoids the redundant spline
  // lookups; it is invalidated in set_physiology when conductance/soil change.
  bool   transpiration_cached_ = false;
  double transpiration_cache_psi_stem_ = 0.0;
  double transpiration_cache_psi_upstream_ = 0.0;
  double transpiration_cache_value_ = 0.0;

  // Cache for the temperature/O2-dependent photosynthesis parameters set in
  // set_physiology (vcmax_, jmax_, gamma_, ko_, kc_, R_d_, km_). They are pure
  // functions of (leaf_temp_, atm_o2_kpa_) and constants, so when the key is
  // unchanged the Arrhenius transcendentals are skipped and the members reused
  // (bit-identical: same inputs -> same outputs). In the current driver
  // leaf_temp_/atm_o2_kpa_ are constant across the run, so this fires once.
  // NOTE: electron_transport_ is deliberately NOT cached here -- it also depends
  // on the per-call PPFD_ and is recomputed every call.
  bool   photo_temp_cached_ = false;
  double photo_temp_cache_leaf_temp_ = 0.0;
  double photo_temp_cache_atm_o2_kpa_ = 0.0;
  std::vector<double> f_r;
  // TODO: move into environment?

  // TODO: atm_vpd - now set in set_physiology although ideally should be moved to enviroment
  double atm_vpd = 2.0; //kPa
  double ca = 40.0; // Pa
  double atm_kpa = 101.3; //kPa
  //partial pressure o2 (kPa)
  double atm_o2_kpa = 21;
  //leaf temperature (deg C)
  double leaf_temp = 25;
  // default soil-moisture content used by the Medlyn beta_ stress factor
  // (matches the fixed values develop's TF24 caller passed to set_physiology)
  double theta_w = 0.2;  //m^3 m^-3
  double theta_fc = 0.5; //m^3 m^-3
  double theta = 0.3;    //m^3 m^-3
  // density of water


  // this might end up hard-coded
  void initialize_integrator(int integration_rule = 21,
                             double integration_tol = 1e-3) {

    integrator = quadrature::QAG(integration_rule,
                                 1, // fixed integration
                                 integration_tol, integration_tol);
  }
  
  // set-up functions
  void set_physiology(double area_leaf, const std::vector<double>& mass_root_prop, double rho, double a_bio, double PPFD, const std::vector<double>& psi_soil, const std::vector<double>& soil_depth, double leaf_specific_conductance_max, double atm_vpd, double ca, double sapwood_volume_per_leaf_area, double leaf_temp, double atm_o2_kpa, double atm_kpa);
  void setup_transpiration(double resolution);
  void setup_root_vulnerability(double resolution);
  // Shared builder for the knot grid {0, step, .., <= psi_max} and the
  // cumulative vulnerability integral G(m) = int_0^m exp(-(s/b)^c) ds, seeded
  // from its gamma closed form. Used by both setup_* functions (see #468).
  void build_cumulative_vulnerability_integral(double b, double c,
                                               double resolution,
                                               std::vector<double>& x,
                                               std::vector<double>& y_integral);
  void setup_clean_leaf();

  // Medlyn stomatal-conductance model (from develop #450); R-callable, standalone.
  double medlyn_model_gs(double assim_colimited_);
  double medlyn_stom_cond_minus_coupled_stom_cond(double x);
  void solve_medlyn_ci_numerical();
  void solve_medlyn_ci_analytical();
  // std::vector<double> root_collar_psi(std::vector<double> soil_moist_);

  void E_from_Soil_to_Root_Collar(double P_x_r, const std::vector<double>& psi_soil);
  void find_root_collar_psi();
  // Shared setup for the root-collar solve: builds the soil-side caches, handles
  // every feasibility early-exit (shutdown / assim<0 / collapsed interval) by
  // setting the final operating point itself, and otherwise returns the feasible
  // collar-potential interval [bound_a, bound_b] (positive magnitudes). Returns
  // false when the operating point is already fully determined (caller is done),
  // true when there is a real interval to choose a collar potential within.
  bool prepare_collar_solve(double& bound_a, double& bound_b);
  // Evaluate the leaf at a *given* root-collar potential (positive magnitude)
  // rather than optimising it: reuses prepare_collar_solve, clamps the target to
  // the feasible interval, and evaluates there (no golden-section search). Leaves
  // exactly the same outputs as find_root_collar_psi and returns profit_. Used by
  // TF24f's gradient-ascent acclimation (#525).
  double evaluate_root_collar_psi(double target_opt_root_psi);
  // Evaluate profit at a given root-collar potential (positive magnitude)
  // *assuming prepare_collar_solve has already run this step* (soil-side caches
  // built, feasible interval [bound_a, bound_b] known). Clamps the target into
  // the interval and sets the operating point (opt_psi_stem_, root_collar_psi_,
  // profit_), returning profit_. This is the post-prepare body of
  // evaluate_root_collar_psi, factored out so the centred finite-difference leaf
  // solve can share one prepare_collar_solve across its three profit evals (#530).
  double profit_at_collar_psi(double target_opt_root_psi,
                              double bound_a, double bound_b);
  // Exact d(profit)/d(opt_root_psi) at a given root-collar potential (positive
  // magnitude), for TF24f's acclimation tracking (#525/#527). Combines
  // forward-mode AD for the analytic photosynthesis/cost algebra, the
  // implicit-function theorem at the psi_stem_to_ci root-find, and analytic
  // spline derivatives (Interpolator::deriv) for the smooth transport. Replaces
  // the noisy finite-difference gradient. Assumes prepare_collar_solve setup has
  // run (psi_soil_inverted_ etc.), as evaluate_root_collar_psi does.
  double dprofit_droot_collar_psi(double opt_root_psi);
  // Analytic d(E_up_)/d(collar potential) for the soil->root-collar uptake
  // (kg H2O m^-2 s^-1 per MPa of signed collar potential P_x_r), mirroring the
  // general branch of E_from_Soil_to_Root_Collar layer by layer. The integral's
  // derivative collapses to +/- root_vuln_integral_from_psi.deriv (the analytic
  // slope of the same pre-integrated vulnerability curve used for the value, so
  // it stays consistent even where that spline extrapolates), so no finite
  // difference is needed.
  // Returns NaN when any layer sits on a branch kink (P_x_r == psi_soil[i], the
  // gravity-balance point, or P_x_r == 0); the caller (dprofit_droot_collar_psi)
  // then falls back to a central difference. Used only on the TF24f acclimation
  // gradient path, not the base TF24 value path.
  double dE_from_soil_dpsi_collar(double P_x_r, const std::vector<double>& psi_soil);
  // Shut-down operating point used by the find_root_collar_psi early-exits: stem
  // held at psi_crit (no transpiration), paying only respiration + hydraulic
  // cost. Only root_collar_psi_ differs between the cases, so it is the argument.
  void set_shutdown_state(double root_collar);
  double find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit);
  double find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil);
  double E_column(double x, const std::vector<double>& psi_soil, double psi_leaf);
  double E_column_zero(double x, const std::vector<double>& psi_soil);
  
  double arrh_curve(double Ea, double ref_value, double leaf_temp) const;
  double peak_arrh_curve(double Ea, double ref_value, double leaf_temp, double H_d, double d_S) const;

  // --- Penman-Monteith leaf energy balance (minimal core; #523) ---------------
  // Recompute the temperature-dependent photosynthetic parameters (vcmax_,
  // jmax_, gamma_, ko_, kc_, R_d_, km_, electron_transport_) at a given leaf
  // temperature. Extracted verbatim from the inline block in set_physiology so
  // the non-PM path is bit-identical; on the PM path it is called per
  // operating point with the energy-balance Tleaf (defeating the photo_temp
  // cache, as intended).
  void update_temperature_dependent_params(double leaf_temp);
  // Saturation vapour pressure es(T) (kPa) and its slope Delta(T) (kPa K^-1),
  // Tetens formula. Not wired into the minimal-cut solve (prescribed atm_vpd is
  // kept); provided and unit-tested for the leaf-to-air VPD in the full cut.
  double saturation_vapour_pressure(double temp) const;
  double saturation_vapour_pressure_slope(double temp) const;
  // Explicit leaf energy balance: Tleaf = Tair + (Rn - lambda*E) * ra / (rho*cp).
  // E is the hydraulically-pinned transpiration (kg H2O m^-2 s^-1); no PM
  // inversion and no A->E feedback, so this is a single algebraic forward pass.
  double leaf_temp_from_E(double E) const;

  // transpiration functions

  // proportion of conductivity in xylem at a given water potential (return: unitless)
  double proportion_of_conductivity(double psi) const;

  // supply-side transpiration for a given water potential gradient between leaves and soil, 
  // references setup_transpiraiton for values (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration(double psi_stem, double psi_upstream);
  // supply-side transpiration for a given water potential gradient between leaves and soil, integrated internally (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration_full_integration(double psi_stem, double psi_upstream);                    
  // stomatal conductance rate of c02 (return: mol CO2 m^-2 s^-1)
  double stom_cond_CO2(double psi_stem, double psi_upstream); // define as a constant
  // converts transpiration in kg h20 s^-1 m^-2 LA to psi_stem (return: -MPa)
  double transpiration_to_psi_stem(double transpiration_, double psi_upstream);
  
  // assimilation functions

  //
  double assim_rubisco_limited(double ci_);
  double electron_transport();
  double assim_electron_limited(double ci_);
  double assim_colimited(double ci_);
  double assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream);
  double psi_stem_to_ci(double psi_stem, double psi_upstream);
  void set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream);


// leaf economics functions
  double hydraulic_cost_Sperry(double psi_stem, double psi_upstream);
  double hydraulic_cost_TF(double psi_stem);

  double profit_psi_stem_Sperry(double psi_stem, double psi_upstream);
  double profit_psi_stem_TF(double psi_stem, double psi_upstream);

// optimiser functions
  void optimise_psi_stem_Sperry();
  void optimise_psi_stem_TF();

};
} // namespace plant
#endif
