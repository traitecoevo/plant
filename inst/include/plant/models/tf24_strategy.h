// Built from  inst/include/plant/models/ff16_strategy.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_STRATEGY_H_
#define PLANT_PLANT_TF24_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/qag.h>
#include <plant/canopy_shape.h>
#include <plant/leaf_model.h>

namespace plant {

class TF24_Strategy: public Strategy<TF24_Environment> {
public:
  typedef std::shared_ptr<TF24_Strategy> ptr;
  TF24_Strategy();

  double compute_average_light_environment(double z, double height,
                                           double height_inverse,
                                           const TF24_Environment &environment);

  // calculate the amount of water transpired relativised by leaf area index.

  double evapotranspiration_dt(double area_leaf_);


  // Overrides ----------------------------------------------

  // Fixed aux slots for the hot ODE rate path, used instead of
  // aux_index.at("...") string-map lookups (see #466). These MUST stay in sync
  // with the order of aux_names() below. If you add/reorder an aux name,
  // update these constants.
  static constexpr int COMPETITION_EFFECT_AUX_INDEX = 0;
  static constexpr int HEIGHT_INVERSE_AUX_INDEX = 1;
  static constexpr int NET_MASS_PRODUCTION_DT_AUX_INDEX = 2;
  static constexpr int AREA_SAPWOOD_AUX_INDEX = 3;

  // update this when the length of state_names changes
  static size_t state_size () { return 5; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
      "area_heartwood",
      "mass_heartwood"
      });
  }

  std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "height_inverse",
      "net_mass_production_dt"
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
  // Inline (header) so it folds into the hot competition path reached from
  // templated Individual<TF24> code (no LTO build).
  double area_leaf(double height) const {
    return std::pow(height / a_l1, 1.0 / a_l2);
  }

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

  double mass_above_ground(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_root) const;

  void compute_rates(const TF24_Environment& environment,
                Internals& vars);

  // Inline (header): per state-set / ODE-state update from templated
  // Individual<TF24> code, avoids a cross-TU call (no LTO build).
  void update_dependent_aux(const int index, Internals& vars) {
    if (index == HEIGHT_INDEX) {
      double height = vars.state(HEIGHT_INDEX);
      vars.set_aux(COMPETITION_EFFECT_AUX_INDEX, area_leaf(height));
      vars.set_aux(HEIGHT_INVERSE_AUX_INDEX, 1.0 / height);
    }
  }

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  double assimilation(const TF24_Environment& environment, double height,
                      double area_leaf, double height_inverse);
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

  // Two overloads: the 3-arg form is the public, virtual entry point used by
  // callers that only have `height` (e.g. establishment_probability); it just
  // computes height_inverse = 1/height and delegates. The 4-arg form is the
  // real (non-virtual) implementation, taking a pre-computed height_inverse so
  // the hot path (compute_rates) can skip the redundant division.
  virtual double net_mass_production_dt(const TF24_Environment& environment,
                                double height, double area_leaf_);
  double net_mass_production_dt(const TF24_Environment& environment,
                                double height, double area_leaf_,
                                double height_inverse);
  // Worker overload that also reports the sapwood intermediates so callers
  // (compute_rates) can reuse them instead of recomputing. Bit-identical: the
  // out-refs receive exactly the values the body already computed.
  double net_mass_production_dt(const TF24_Environment& environment,
                                double height, double area_leaf_,
                                double height_inverse,
                                double& area_sapwood_, double& mass_sapwood_);

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


  double mortality_dt(double productivity_area, double cumulative_mortality) const;
  double mortality_growth_independent_dt()const ;
  double mortality_growth_dependent_dt(double productivity_area) const;
  // [eqn 20] Survival of seedlings during establishment
  double establishment_probability(const TF24_Environment& environment);

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  // Inline (header) so the per-node hot competition path called from
  // Individual<TF24>::compute_competition inlines these helpers instead of
  // paying a cross-TU call each iteration (no LTO build).
  double compute_competition(double z, double height) const {
    return compute_competition(z, area_leaf(height), 1.0 / height);
  }
  double compute_competition(double z, double area_leaf_,
                             double height_inverse) const {
    return compute_competition_by_ratio(z * height_inverse, area_leaf_);
  }
  double compute_competition_by_ratio(double z_over_height,
                                      double area_leaf_) const {
    return k_I * area_leaf_ * canopy_shape.Q(z_over_height);
  }

  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  double Qp(double x, double height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // Set constants within TF24_Strategy
  void prepare_strategy();

  // * Core traits
  double lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  double rho       = 608.0;      // Wood density [kg/m3]
  double hmat      = 16.5958691; // Height at maturation [m]
  double omega     = 3.8e-5;     // Seed mass [kg]
  // * Individual allometry
  // Canopy shape parameters
  double eta       = 12.0; // [dimensionless]
  double eta_c     = NA_REAL; // [dimensionless]
  CanopyShape canopy_shape;
  // Sapwood area per leaf area
  // Ratio sapwood area area to leaf area
  double theta     = 1.0/4669; // [dimensionless]
  // Height - leaf mass scaling
  double a_l1        = 5.44; // height with 1m2 leaf [m]
  double a_l2        = 0.306; // dimensionless scaling of height with leaf area
  // Root mass per leaf area
  double a_r1        = 0.07;  //[kg / m]
  // Ratio of bark area : sapwood area
  double a_b1         = 0.17; // [dimensionless]

  // * Production
  // Ratio of leaf dark respiration to leaf mass [mol CO2 / yr  / kg]
  // =  [mol CO2 / m2 / yr]  |  (39.27 = 2100 * 0.00187)  | narea * photosynthesis_per_nitrogen
  //    / [kg(leaf) / m2 ]   |    / (0.1978791)           | lma
  // Hard coded in value of lma here so that this value doesn't change
  // if that trait changes above.
  double r_l    = 39.27 / 0.1978791;
  // Root respiration per mass [mol CO2 / yr / kg]
  double r_r    = 217.0;
  // Sapwood respiration per stem mass  [mol CO2 / yr / kg]
  // = respiration per volume [mol CO2 / m3 / yr]
  // /  wood density [kg/m3]
  double r_s    = 4012.0 / 608.0;
  // Bark respiration per stem mass
  // assumed to be twice rate of sapwood
  // (NOTE that there is a re-parametrisation here relative to the paper
  // -- r_b is defined (new) as 2*r_s, whereas the paper assumes a
  // fixed multiplication by 2)
  double r_b    = 2.0 * r_s;
  // Carbon conversion parameter
  double a_y    = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  double a_bio  = 2.45e-2;
  // Leaf turnover [/yr]
  double k_l    =  0.4565855;
  // Bark turnover [/yr]
  double k_b    = 0.2;
  // Sapwood turnover [/yr]
  double k_s           = 0.2;
  // Root turnover [/yr]
  double k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  double a_p1   = 151.177775377968; // [mol CO2 / yr / m2]
  double a_p2   = 0.204716166503633; // [dimensionless]

  // * Seed production
  // Accessory cost of reproduction
  double a_f3  = 3.0 *  3.8e-5; // [kg per seed]
  // Maximum allocation to reproduction
  double a_f1   = 1.0; //[dimensionless]
  // Size range across which individuals mature
  double a_f2   = 50; // [dimensionless]

  // * Mortality parameters
  // Probability of survival during dispersal
  double S_D   = 0.25; // [dimensionless]
  // Parameter for seedling survival
  double a_d0    = 0.1; //[kg / yr / m2]
  // Baseline for intrinsic mortality
  double d_I    = 0.01; // [ / yr]
  // Baseline rate for growth-related mortality
  double a_dG1    = 5.5; // [ / yr]
  // Risk coefficient for dry mass production (per area)
  double a_dG2    = 20.0;// [yr m2 / kg ]

  // Germination
  double recruitment_decay = 0.0;
  
  // * Light capture parameters
  double k_I = 0.5;

  // Height and leaf area of a (germinated) seed
  double height_0  = NA_REAL;
  double height_0_inverse = NA_REAL;
  double area_leaf_0;


  Leaf leaf;

  // leaf traits - default values Eucalyptus saligna
  double vcmax_25 = 96;
  double p_50 = 1.85;
  double K_s = 1;
  double c = log(log(1-0.5)/log(1-0.88))/(log(p_50) - log(5.16));
  double b = p_50 / std::pow(-log(1 - 50.0 / 100.0), 1 / c);
  double psi_crit = b*std::pow(log(1/0.05),1/c); // derived from b and c
  double beta1 = 20000;
  double beta2 = 1.5;
  double jmax_25 = vcmax_25*1.64;
  double hk_s = 4;
  double a = 0.30; // effective quantum yield of electron transport  (mol photon mol ^-1 electron)  Sabot et al. 2020
  double curv_fact_elec_trans = 0.7; 
  double curv_fact_colim = 0.99; 
  double var_sapwood_volume_cost = 1; 
  double newton_tol_abs = 0.001;
  double GSS_tol_abs = 1e-7;
  double vulnerability_curve_ncontrol = 100;
  double ci_abs_tol = 1e-6;
  double ci_niter = 1000;
  double g0 = 0.022; //g0 parameter in the medlyn model umol m^-2 s^-1
  double g1 = 2.57; //g1 parameter in the medlyn model umol kPa^0.5
  double g1_TF24 = 46.32995;

  //nitrogen allocation traits (parameterised from Austraits 4.1.0)
  double nmass_l = 13e-3; // kg N kg^-1 mass
  double nmass_s = 1.98e-3; // kg N kg^-1 mass
  double nmass_b = 3.40e-3; // kg N kg^-1 mass
  double nmass_r = 3.35e-3; // kg N kg^-1 mass
  double dmass_dN = 0; //change in mass per change in kg kg^-1 N

  std::string name;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;
};

TF24_Strategy::ptr make_strategy_ptr(TF24_Strategy s);

}

#endif
