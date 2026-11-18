// -*-c++-*-
#ifndef PLANT_PLANT_ES20_MODEL_H_
#define PLANT_PLANT_ES20_MODEL_H_

// ES20: FF16 with a non-structural carbohydrate (NSC) storage pool, after
// Stefaniak, Tissue, Falster & Medlyn (egusphere-2026-1474; called plantNSC
// there). Net carbon uptake enters storage first; growth and reproduction draw
// on storage. Two traits set the allocation strategy: the switch time `t_s`,
// the time of year after which a plant stops drawing on storage to grow, and
// the utilisation rate `a_s`. Each year has a stress period, from a stress
// onset to the end of the year, during which photosynthesis, establishment and
// productivity-driven mortality are switched off and respiration is reduced.
//
// Ported from Elisa Stefaniak's ES20_Strategy (foxeswithdata/plant @ 541ba5ed).
// The equations are hers, including two that look like slips and are kept so
// that the port reproduces her runs (see "Kept as written" in compute_rates).
// What changed is the plumbing:
//  - The stress schedule is the "stress_onset" extrinsic driver on the
//    environment (one value per year, the fraction of the year at which stress
//    starts), so ES20 runs in FF16's light environment rather than its own.
//    With no "stress_onset" driver there is no stress at all.
//  - Parameters live in ES20_Pars, with the values of Table S2 of the paper.
//  - Crown assimilation uses plant's fixed crown quadrature, not adaptive QAG.
//  - Her constants 40 (calibration rate), 50 (stress transition) and 0.5
//    (respiration kept during stress) are the parameters c_k, k_stress, r_min.

#include <plant/models/model.h>
#include <plant/models/ff16_environment.h>
#include <plant/models/canopy_shape.h>
#include <plant/numerics/qk.h>
#include <plant/numerics/with_slope.h>
#include <odelia/root_find.hpp>
#include <cmath>
#include <limits>

namespace plant {

// Biological parameters. Defaults are Table S2 of the paper; where a 2023
// ES20_Strategy default differed, it is given in the comment.
struct ES20_Pars {
  // * Core traits
  double lma   = 0.1242302;  // Leaf mass per area [kg / m2]
  double rho   = 505.0;      // Wood density [kg / m3]
  double hmat  = 30.0;       // Height at maturation [m]
  double omega = 0.1;        // Seed mass [kg]
  // * Allometry. Leaf area is exp(a_l1 + a_l2 log(height)), not FF16's form.
  double eta   = 5.0;          // Crown shape [-] (2023 code: 12)
  double theta = 0.000150887;  // Sapwood area per leaf area [-]
  double a_l1  = -0.56936;     // Leaf area-height intercept (log scale)
  double a_l2  = 1.48185;      // Leaf area-height exponent [-]
  double a_r1  = 0.027691;     // Root mass per leaf area [kg / m2]
  double a_b1  = 0.3768273;    // Bark area per sapwood area [-]
  // * Production
  double r_l   = 175.8027;     // Leaf respiration per mass [mol / yr / kg] (2023 code: 39.27 / 0.1978791)
  double r_r   = 217.0;        // Root respiration per mass [mol / yr / kg]
  double r_s   = 2.412504;     // Sapwood respiration per mass [mol / yr / kg]
  double r_b   = 2.0 * 2.412504; // Bark respiration per mass [mol / yr / kg]
  double r_min = 0.5;          // Fraction of respiration kept under full stress [-]
  double a_y   = 0.7;          // Yield [-]
  double a_bio = 2.45e-2;      // Dry mass per mol carbon [kg / mol]
  double k_l   = 1.012133;     // Leaf turnover [/ yr] (2023 code: 0.4565855)
  double k_b   = 0.06;         // Bark turnover [/ yr]
  double k_s   = 0.06;         // Sapwood turnover [/ yr]
  double k_r   = 1.0;          // Root turnover [/ yr]
  double a_p1  = 82.34303;     // Leaf photosynthesis per area [mol / yr / m2] (2023 code: 151.18)
  double a_p2  = 0.1072532;    // Saturation of leaf photosynthesis [-] (2023 code: 0.2047)
  // * Reproduction
  double a_f1  = 0.05;         // Maximum allocation to reproduction [-] (2023 code: 0.10)
  double a_f2  = 50.0;         // Steepness of maturation around hmat [-]
  double a_f3  = 3.0 * 3.8e-5; // Accessory cost per seed [kg]
  // * Mortality and establishment
  double S_D   = 0.25;   // Survival through dispersal [-]
  double a_d0  = 0.1;    // Establishment parameter [kg / yr / m2]
  double d_I   = 0.01;   // Growth-independent mortality [/ yr]
  double a_dP1 = 0.25;   // Productivity-dependent mortality, baseline [/ yr] (2023 code a_dG1: 3)
  double a_dP2 = 20.0;   // Productivity-dependent mortality, risk [yr m2 / kg] (2023 code a_dG2)
  double a_dS1 = 3.0;    // Storage-dependent mortality, baseline [/ yr]
  double a_dS2 = 300.0;  // Storage-dependent mortality, risk [-]
  // * Storage
  double t_s   = 0.4109589;  // Switch time, fraction of the year [-] (paper: 0.25 or 0.5)
  double a_s   = 0.06 * 365; // Storage utilisation rate [/ yr] (paper: 0.1 or 0.3 / day)
  double b_s1  = 0.1;        // Initial storage per unit live mass [-]
  double b_s2  = 120.0;      // Steepness of the growth switch [/ yr]
  double c_k   = 40.0;       // Allometric calibration rate (Eq. 3) [per unit departure]
  // * Environment
  double k_stress = 50.0;    // Steepness of the stress transitions [/ yr]
  double k_I   = 0.5;        // Light extinction coefficient [-]
};

class ES20_Model : public Model<FF16_Environment<double>> {
public:
  using value_type = double;
  typedef FF16_Environment<double> environment_type;
  typedef std::shared_ptr<ES20_Model> ptr;
  ES20_Model();

  // Scientific version. Bump only when equations or default parameters change
  // the simulation output for identical inputs (see agents.md §7).
  static constexpr int scientific_version = 1;

  // Growth reads the storage pool, so it is not a function of size alone and
  // the height coordinate's compression term does not hold.
  static constexpr bool density_in_birth_date_default = true;

  // State slots. The order is Elisa's, and the first five are FF16's. MUST stay
  // in sync with state_names().
  static constexpr int AREA_HEARTWOOD_INDEX = 3;
  static constexpr int MASS_HEARTWOOD_INDEX = 4;
  static constexpr int MASS_STORAGE_INDEX = 5;
  static constexpr int AREA_LEAF_INDEX = 6;
  static constexpr int MASS_LEAF_INDEX = 7;
  static constexpr int AREA_SAPWOOD_INDEX = 8;
  static constexpr int MASS_SAPWOOD_INDEX = 9;
  static constexpr int AREA_BARK_INDEX = 10;
  static constexpr int MASS_BARK_INDEX = 11;
  static constexpr int MASS_ROOT_INDEX = 12;
  static constexpr int AREA_STEM_INDEX = 13;
  static constexpr int DIAMETER_STEM_INDEX = 14;

  static size_t state_size() { return 15; }
  static std::vector<std::string> state_names() {
    return {"height", "mortality", "fecundity", "area_heartwood",
            "mass_heartwood", "mass_storage", "area_leaf", "mass_leaf",
            "area_sapwood", "mass_sapwood", "area_bark", "mass_bark",
            "mass_root", "area_stem", "diameter_stem"};
  }

  // Aux slots. MUST stay in sync with aux_names().
  enum Aux : int {
    COMPETITION_EFFECT_AUX_INDEX = 0, HEIGHT_INVERSE_AUX_INDEX,
    NET_MASS_PRODUCTION_DT_AUX_INDEX, DBIOMASS_DT_AUX_INDEX,
    RESPIRATION_DT_AUX_INDEX, RESPIRATION_LEAF_DT_AUX_INDEX,
    RESPIRATION_SAPWOOD_DT_AUX_INDEX, RESPIRATION_ROOT_DT_AUX_INDEX,
    RESPIRATION_BARK_DT_AUX_INDEX, DAREA_LEAF_DMASS_LIVE_AUX_INDEX,
    DAREA_LEAF_DT_AUX_INDEX, DMASS_LEAF_DAREA_LEAF_AUX_INDEX,
    DMASS_SAPWOOD_DAREA_LEAF_AUX_INDEX, DMASS_BARK_DAREA_LEAF_AUX_INDEX,
    DMASS_ROOT_DAREA_LEAF_AUX_INDEX, DMASS_LEAF_DT_AUX_INDEX,
    DMASS_SAPWOOD_DT_AUX_INDEX, DMASS_BARK_DT_AUX_INDEX,
    DMASS_ROOT_DT_AUX_INDEX, DHEIGHT_DAREA_LEAF_AUX_INDEX,
    HEIGHT_DIFFERENCE_AUX_INDEX, HEIGHT_ADJUSTMENT_AUX_INDEX,
    MASS_SAPWOOD_DIFFERENCE_AUX_INDEX, MASS_SAPWOOD_ADJUSTMENT_AUX_INDEX,
    FECUNDITY_DT_ABS_AUX_INDEX, STORAGE_PORTION_AUX_INDEX, DEAD_FLAG_AUX_INDEX,
    MORTALITY_STORAGE_DEPENDENT_DT_AUX_INDEX,
    MORTALITY_PRODUCTIVITY_DEPENDENT_DT_AUX_INDEX, PRODUCTIVITY_AREA_AUX_INDEX,
    STRESS_AUX_INDEX, AUX_COUNT
  };
  size_t aux_size() { return AUX_COUNT; }
  std::vector<std::string> aux_names() {
    return {"competition_effect", "height_inverse", "net_mass_production_dt",
            "dbiomass_dt", "respiration_dt", "respiration_leaf_dt",
            "respiration_sapwood_dt", "respiration_root_dt",
            "respiration_bark_dt", "darea_leaf_dmass_live", "darea_leaf_dt",
            "dmass_leaf_darea_leaf", "dmass_sapwood_darea_leaf",
            "dmass_bark_darea_leaf", "dmass_root_darea_leaf", "dmass_leaf_dt",
            "dmass_sapwood_dt", "dmass_bark_dt", "dmass_root_dt",
            "dheight_darea_leaf", "height_difference", "height_adjustment",
            "mass_sapwood_difference", "mass_sapwood_adjustment",
            "fecundity_dt_abs", "storage_portion", "dead_flag",
            "mortality_storage_dependent_dt",
            "mortality_productivity_dependent_dt", "productivity_area",
            "stress"};
  }

  void refresh_indices();

  // * Environment
  // Stress factor E_s in [0, 1]: ~1 before the year's stress onset, ~0 after
  // it, with logistic transitions of steepness k_stress (Elisa's getStress()).
  double stress(const environment_type& environment) const;
  // The growth switch: ~1 before t_s in each year, ~0 after (her dbiomass_dt).
  double growth_switch(double time) const;

  // * Allometry (Elisa's equations)
  double area_leaf(double height) const {
    return std::exp(pars.a_l1 + pars.a_l2 * std::log(height));
  }
  double mass_leaf(double area_leaf) const { return area_leaf * pars.lma; }
  double area_sapwood(double area_leaf) const { return area_leaf * pars.theta; }
  double mass_sapwood(double area_sapwood, double height) const {
    return area_sapwood * height * eta_c * pars.rho;
  }
  double area_bark(double area_leaf) const {
    return pars.a_b1 * area_leaf * pars.theta;
  }
  double mass_bark(double area_bark, double height) const {
    return area_bark * height * eta_c * pars.rho;
  }
  double area_stem(double area_bark, double area_sapwood,
                   double area_heartwood) const {
    return area_bark + area_sapwood + area_heartwood;
  }
  double diameter_stem(double area_stem) const {
    return std::sqrt(4 * area_stem / M_PI);
  }
  double mass_root(double area_leaf) const { return pars.a_r1 * area_leaf; }
  double mass_live(double mass_leaf, double mass_bark, double mass_sapwood,
                   double mass_root) const {
    return mass_leaf + mass_sapwood + mass_bark + mass_root;
  }
  double mass_live_given_height(double height) const;
  double height_given_mass_leaf(double mass_leaf) const;

  // Allometric calibration, Eqs. 2-3: a factor in (0, 2) that speeds the growth
  // of a component that is below its allometric size and slows one above it.
  double calibration(double difference) const {
    return 2.0 / (1.0 + std::exp(-pars.c_k * difference));
  }
  double dheight_darea_leaf(double area_leaf, double height) const;
  double dmass_sapwood_darea_leaf(double area_leaf, double height,
                                  double mass_sapwood) const;
  double dmass_bark_darea_leaf(double area_leaf, double height,
                               double mass_sapwood) const {
    return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf, height, mass_sapwood);
  }
  double darea_leaf_dmass_live(double area_leaf, double height,
                               double mass_sapwood) const {
    return 1.0 / (pars.lma +
                  dmass_sapwood_darea_leaf(area_leaf, height, mass_sapwood) +
                  dmass_bark_darea_leaf(area_leaf, height, mass_sapwood) +
                  pars.a_r1);
  }

  // * Production
  double assimilation(const environment_type& environment, double height,
                      double area_leaf, double stress) const;
  double respiration(double mass_leaf, double mass_sapwood, double mass_bark,
                     double mass_root, double stress) const;
  double net_mass_production_dt(const environment_type& environment,
                                double height, double area_leaf,
                                double mass_leaf, double mass_sapwood,
                                double mass_bark, double mass_root) const;
  double net_mass_production_dt(const environment_type& environment,
                                const Internals<double>& vars) {
    return net_mass_production_dt(
        environment, vars.state(HEIGHT_INDEX), vars.state(AREA_LEAF_INDEX),
        vars.state(MASS_LEAF_INDEX), vars.state(MASS_SAPWOOD_INDEX),
        vars.state(MASS_BARK_INDEX), vars.state(MASS_ROOT_INDEX));
  }

  // * Storage, reproduction and mortality
  // Utilisation U (Eq. 1): what the plant draws from storage to grow.
  double dbiomass_dt(double time, double mass_storage) const;
  double fecundity_dt(double mass_storage, double height) const;
  double mortality_dt(double storage_portion, double productivity_area,
                      double stress, double cumulative_mortality) const;
  double mortality_storage_dependent_dt(double storage_portion) const {
    return pars.a_dS1 * std::exp(-pars.a_dS2 * storage_portion);
  }
  double mortality_productivity_dependent_dt(double productivity_area,
                                             double stress) const {
    return stress * pars.a_dP1 * std::exp(-pars.a_dP2 * productivity_area);
  }

  void compute_rates(const environment_type& environment, Internals<double>& vars);

  // * Birth
  double initial_height() const { return height_0; }
  // Writes the seedling's tissues and storage pool (Eq. 8). The core calls it
  // for every node and stochastic arrival; for a single plant, ES20_Individual()
  // in R calls it.
  void set_initial_states(const environment_type& environment,
                          Internals<double>& vars);
  double establishment_probability(const environment_type& environment);
  double establishment_probability(const environment_type& environment,
                                   const Internals<double>& vars) {
    return establishment_probability_given(
        vars.aux(NET_MASS_PRODUCTION_DT_AUX_INDEX));
  }
  double establishment_probability_given(double net_mass_production_dt) const;

  // * Competition: the leaf area a plant actually holds, a state, not the
  // allometric leaf area of its height. Inline: this is the per-node hot path.
  void update_dependent_aux(const int index, Internals<double>& vars) {
    if (index == HEIGHT_INDEX) {
      vars.set_aux(HEIGHT_INVERSE_AUX_INDEX, 1.0 / vars.state(HEIGHT_INDEX));
    } else if (index == AREA_LEAF_INDEX) {
      vars.set_aux(COMPETITION_EFFECT_AUX_INDEX, vars.state(AREA_LEAF_INDEX));
    }
  }
  double compute_competition(double z, const Internals<double>& vars) const {
    return pars.k_I * vars.aux(COMPETITION_EFFECT_AUX_INDEX) *
           canopy_shape.leaf_area_above(z * vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }
  with_slope<double>
  compute_competition_and_slope(double z, const Internals<double>& vars) const {
    const double scale = pars.k_I * vars.aux(COMPETITION_EFFECT_AUX_INDEX);
    const double height_inverse = vars.aux(HEIGHT_INVERSE_AUX_INDEX);
    const std::pair<double, double> Qq =
        canopy_shape.Q_and_q(z * height_inverse, z, height_inverse);
    return {scale * Qq.first, -(scale * Qq.second)};
  }

  double height_seed() const;
  void prepare_strategy();

  ES20_Pars pars;

  // Derived in prepare_strategy() (not user-set)
  double eta_c = NA_REAL;
  CanopyShape<double> canopy_shape;
  quadrature::QK function_integrator;
  crown_rule crown;
  double height_0 = NA_REAL;
  double area_leaf_0 = NA_REAL, mass_leaf_0 = NA_REAL,
         area_sapwood_0 = NA_REAL, mass_sapwood_0 = NA_REAL,
         area_bark_0 = NA_REAL, mass_bark_0 = NA_REAL, mass_root_0 = NA_REAL,
         area_stem_0 = NA_REAL, diameter_stem_0 = NA_REAL,
         mass_storage_0 = NA_REAL;
};

// Definitions ----------------------------------------------------------------
// Header-only, as K93_Model is.

inline ES20_Model::ES20_Model() {
  collect_all_auxiliary = false;
  refresh_indices();
  name = "ES20";
}

inline void ES20_Model::refresh_indices() {
  state_index = std::map<std::string, int>();
  aux_index = std::map<std::string, int>();
  const std::vector<std::string> s = state_names(), a = aux_names();
  for (size_t i = 0; i < s.size(); i++) state_index[s[i]] = i;
  for (size_t i = 0; i < a.size(); i++) aux_index[a[i]] = i;
  check_state_layout(state_index, state_size(), "ES20");
  if (aux_index.size() != static_cast<size_t>(AUX_COUNT)) {
    util::stop("ES20 names " + util::to_string(static_cast<int>(aux_index.size())) +
               " aux slots against the enum's " +
               util::to_string(static_cast<int>(AUX_COUNT)));
  }
}

// Elisa's piecewise logistic: each branch is a logistic in time centred on one
// transition, and the branches meet halfway between transitions, where the
// logistic is within exp(-k * gap / 2) of 0 or 1.
inline double ES20_Model::stress(const environment_type& environment) const {
  const odelia::drivers::Function* onset_driver =
      environment.extrinsic_drivers.get_function_ptr("stress_onset");
  if (onset_driver == nullptr) {
    return 1.0;
  }
  const double time = environment.time;
  const double yr = std::floor(time);
  const double t_crit = yr + onset_driver->evaluate(yr);
  const double k = pars.k_stress;
  if (time < (yr + t_crit) / 2) {
    return 1 / (1 + std::exp(-k * (time - yr)));
  } else if (time < (t_crit + yr + 1) / 2) {
    return 1 / (1 + std::exp(k * (time - t_crit)));
  }
  return 1 / (1 + std::exp(-k * (time - (yr + 1))));
}

inline double ES20_Model::growth_switch(double time) const {
  const double yr = std::floor(time);
  const double t_s = pars.t_s, b = pars.b_s2;
  if (time < yr + t_s / 2) {
    return 1 / (1 + std::exp(-b * (time - yr)));
  } else if (time < yr + t_s + (1 - t_s) / 2) {
    return 1 / (1 + std::exp(b * (time - (t_s + yr))));
  }
  return 1 / (1 + std::exp(-b * (time - (yr + 1))));
}

inline double ES20_Model::dbiomass_dt(double time, double mass_storage) const {
  return mass_storage > 0 ? pars.a_s * mass_storage * growth_switch(time) : 0.0;
}

inline double ES20_Model::fecundity_dt(double mass_storage, double height) const {
  return mass_storage * pars.a_f1 /
         ((pars.omega + pars.a_f3) *
          (1.0 + std::exp(pars.a_f2 * (1.0 - height / pars.hmat))));
}

inline double ES20_Model::mass_live_given_height(double height) const {
  const double al = area_leaf(height);
  return mass_leaf(al) + mass_bark(area_bark(al), height) +
         mass_sapwood(area_sapwood(al), height) + mass_root(al);
}

inline double ES20_Model::height_given_mass_leaf(double mass_leaf) const {
  return std::exp(1 / pars.a_l2 * std::log(mass_leaf / pars.lma) -
                  pars.a_l1 / pars.a_l2);
}

inline double ES20_Model::dheight_darea_leaf(double area_leaf,
                                             double height) const {
  const double adjustment = calibration(area_leaf - this->area_leaf(height));
  return (std::exp(std::log(area_leaf) / pars.a_l2 - pars.a_l1 / pars.a_l2) /
          (area_leaf * pars.a_l2)) *
         adjustment;
}

inline double ES20_Model::dmass_sapwood_darea_leaf(double area_leaf,
                                                   double height,
                                                   double mass_sapwood) const {
  const double adjustment = calibration(
      this->mass_sapwood(area_sapwood(area_leaf), height) - mass_sapwood);
  return pars.rho * eta_c * pars.theta *
         (height + area_leaf * dheight_darea_leaf(area_leaf, height)) *
         adjustment;
}

// Gross annual assimilation, scaled by the stress factor. Deep-crown form:
// leaf photosynthesis integrated over crown depth.
inline double ES20_Model::assimilation(const environment_type& environment,
                                       double height, double area_leaf,
                                       double stress) const {
  const double a_p1 = pars.a_p1, a_p2 = pars.a_p2;
  return area_leaf * stress *
         environment.light_availability.crown_integral(
             height, crown,
             [a_p1, a_p2](double E) { return a_p1 * E / (E + a_p2); });
}

inline double ES20_Model::respiration(double mass_leaf, double mass_sapwood,
                                      double mass_bark, double mass_root,
                                      double stress) const {
  const double f = pars.r_min + (1 - pars.r_min) * stress;
  return (pars.r_l * mass_leaf + pars.r_b * mass_bark +
          pars.r_s * mass_sapwood + pars.r_r * mass_root) * f;
}

// Net mass production B (dB/dt in the paper). Unlike FF16, turnover is not
// subtracted here: each tissue's turnover is a decay term on its own state.
inline double ES20_Model::net_mass_production_dt(
    const environment_type& environment, double height, double area_leaf,
    double mass_leaf, double mass_sapwood, double mass_bark,
    double mass_root) const {
  const double E_s = stress(environment);
  return pars.a_bio * pars.a_y *
         (assimilation(environment, height, area_leaf, E_s) -
          respiration(mass_leaf, mass_sapwood, mass_bark, mass_root, E_s));
}

inline double ES20_Model::mortality_dt(double storage_portion,
                                       double productivity_area, double stress,
                                       double cumulative_mortality) const {
  if (!mortality_hazard_live(cumulative_mortality)) {
    return 0.0;
  }
  return pars.d_I + mortality_storage_dependent_dt(storage_portion) +
         mortality_productivity_dependent_dt(productivity_area, stress);
}

inline void ES20_Model::compute_rates(const environment_type& environment,
                                      Internals<double>& vars) {
  const double height = vars.state(HEIGHT_INDEX);
  // Her floor on leaf area, applied to the rates only: the 2023 code wrote it
  // back into the state from inside the rate evaluation. Turnover alone takes
  // leaf area to zero only asymptotically, so it does not bind in practice.
  const double area_leaf_ = std::max(vars.state(AREA_LEAF_INDEX), 1e-4);
  const double area_bark_ = vars.state(AREA_BARK_INDEX);
  const double area_heartwood_ = vars.state(AREA_HEARTWOOD_INDEX);
  const double area_sapwood_ = vars.state(AREA_SAPWOOD_INDEX);
  const double mass_leaf_ = vars.state(MASS_LEAF_INDEX);
  const double mass_bark_ = vars.state(MASS_BARK_INDEX);
  const double mass_sapwood_ = vars.state(MASS_SAPWOOD_INDEX);
  const double mass_root_ = vars.state(MASS_ROOT_INDEX);
  const double mass_storage_ = vars.state(MASS_STORAGE_INDEX);

  const double E_s = stress(environment);
  const double respiration_f = pars.a_bio * pars.a_y * (pars.r_min + (1 - pars.r_min) * E_s);
  const double net_mass_production_dt_ =
      pars.a_bio * pars.a_y *
      (assimilation(environment, height, area_leaf_, E_s) -
       respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_, E_s));
  const double dbiomass_dt_ = dbiomass_dt(environment.time, mass_storage_);
  double dmass_storage_dt_ = net_mass_production_dt_ - dbiomass_dt_;
  const double storage_portion_ =
      mass_storage_ / mass_live(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  const double productivity_area_ = net_mass_production_dt_ / area_leaf_;

  vars.set_aux(NET_MASS_PRODUCTION_DT_AUX_INDEX, net_mass_production_dt_);
  vars.set_aux(DBIOMASS_DT_AUX_INDEX, dbiomass_dt_);
  vars.set_aux(STORAGE_PORTION_AUX_INDEX, storage_portion_);
  vars.set_aux(DEAD_FLAG_AUX_INDEX, mass_storage_ < 0 ? 1.0 : 0.0);
  vars.set_aux(RESPIRATION_LEAF_DT_AUX_INDEX, respiration_f * pars.r_l * mass_leaf_);
  vars.set_aux(RESPIRATION_BARK_DT_AUX_INDEX, respiration_f * pars.r_b * mass_bark_);
  vars.set_aux(RESPIRATION_ROOT_DT_AUX_INDEX, respiration_f * pars.r_r * mass_root_);
  vars.set_aux(RESPIRATION_SAPWOOD_DT_AUX_INDEX, respiration_f * pars.r_s * mass_sapwood_);
  vars.set_aux(RESPIRATION_DT_AUX_INDEX,
               respiration_f * (pars.r_l * mass_leaf_ + pars.r_b * mass_bark_ +
                                pars.r_s * mass_sapwood_ + pars.r_r * mass_root_));
  const double height_difference = area_leaf_ - area_leaf(height);
  const double mass_sapwood_difference =
      mass_sapwood(area_sapwood(area_leaf_), height) - mass_sapwood_;
  vars.set_aux(HEIGHT_DIFFERENCE_AUX_INDEX, height_difference);
  vars.set_aux(HEIGHT_ADJUSTMENT_AUX_INDEX, calibration(height_difference));
  vars.set_aux(MASS_SAPWOOD_DIFFERENCE_AUX_INDEX, mass_sapwood_difference);
  vars.set_aux(MASS_SAPWOOD_ADJUSTMENT_AUX_INDEX, calibration(mass_sapwood_difference));

  if (dbiomass_dt_ > 0) {
    const double dheight_darea_leaf_ = dheight_darea_leaf(area_leaf_, height);
    const double darea_leaf_dmass_live_ =
        darea_leaf_dmass_live(area_leaf_, height, mass_sapwood_);
    const double darea_leaf_dt_ = darea_leaf_dmass_live_ * dbiomass_dt_;
    const double dmass_leaf_dt_ = pars.lma * darea_leaf_dt_;

    vars.set_aux(DHEIGHT_DAREA_LEAF_AUX_INDEX, dheight_darea_leaf_);
    vars.set_aux(DAREA_LEAF_DMASS_LIVE_AUX_INDEX, darea_leaf_dmass_live_);
    vars.set_aux(DAREA_LEAF_DT_AUX_INDEX, darea_leaf_dt_);
    vars.set_aux(DMASS_LEAF_DT_AUX_INDEX, dmass_leaf_dt_);
    vars.set_aux(DMASS_LEAF_DAREA_LEAF_AUX_INDEX, pars.lma);

    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf_ * darea_leaf_dt_);
    vars.set_rate(AREA_LEAF_INDEX, darea_leaf_dt_ - pars.k_l * area_leaf_);
    vars.set_rate(MASS_LEAF_INDEX, dmass_leaf_dt_ - pars.k_l * mass_leaf_);

    // Kept as written: heartwood area gains k_s * theta * area_sapwood here,
    // because her area_heartwood_dt() passes sapwood area to area_sapwood(),
    // which multiplies by theta again, while sapwood loses k_s * area_sapwood.
    // The no-growth branch below adds k_s * area_sapwood. Heartwood area reaches
    // only diameter_stem, which nothing else reads.
    const double darea_heartwood_dt_ = pars.k_s * area_sapwood(area_sapwood_);
    vars.set_rate(AREA_HEARTWOOD_INDEX, darea_heartwood_dt_);
    vars.set_rate(MASS_HEARTWOOD_INDEX, pars.k_s * mass_sapwood_);

    const double darea_sapwood_dt_ = darea_leaf_dt_ * pars.theta;
    const double dmass_sapwood_darea_leaf_ =
        dmass_sapwood_darea_leaf(area_leaf_, height, mass_sapwood_);
    const double dmass_sapwood_dt_ = darea_leaf_dt_ * dmass_sapwood_darea_leaf_;
    vars.set_aux(DMASS_SAPWOOD_DAREA_LEAF_AUX_INDEX, dmass_sapwood_darea_leaf_);
    vars.set_aux(DMASS_SAPWOOD_DT_AUX_INDEX, dmass_sapwood_dt_);
    vars.set_rate(AREA_SAPWOOD_INDEX, darea_sapwood_dt_ - pars.k_s * area_sapwood_);
    vars.set_rate(MASS_SAPWOOD_INDEX, dmass_sapwood_dt_ - pars.k_s * mass_sapwood_);

    const double darea_bark_dt_ = pars.a_b1 * darea_leaf_dt_ * pars.theta;
    const double dmass_bark_dt_ = pars.a_b1 * dmass_sapwood_dt_;
    vars.set_aux(DMASS_BARK_DAREA_LEAF_AUX_INDEX, dmass_bark_dt_ / darea_leaf_dt_);
    vars.set_aux(DMASS_BARK_DT_AUX_INDEX, dmass_bark_dt_);
    vars.set_rate(AREA_BARK_INDEX, darea_bark_dt_ - pars.k_b * area_bark_);
    vars.set_rate(MASS_BARK_INDEX, dmass_bark_dt_ - pars.k_b * mass_bark_);

    // Stem area: sapwood, heartwood and bark growth less bark turnover (sapwood
    // turnover moves to heartwood). Diameter grows with the sum of the three
    // component areas, not with the area_stem state.
    const double darea_stem_dt_ = darea_sapwood_dt_ + darea_heartwood_dt_ +
                                  darea_bark_dt_ - pars.k_b * area_bark_;
    const double area_stem_sum = area_stem(area_heartwood_, area_sapwood_, area_bark_);
    vars.set_rate(DIAMETER_STEM_INDEX,
                  std::pow(M_PI * area_stem_sum, -0.5) * darea_stem_dt_);
    vars.set_rate(AREA_STEM_INDEX, darea_stem_dt_);

    const double dmass_root_dt_ = darea_leaf_dt_ * pars.a_r1;
    vars.set_aux(DMASS_ROOT_DAREA_LEAF_AUX_INDEX, dmass_root_dt_ / darea_leaf_dt_);
    vars.set_aux(DMASS_ROOT_DT_AUX_INDEX, dmass_root_dt_);
    vars.set_rate(MASS_ROOT_INDEX, dmass_root_dt_ - pars.k_r * mass_root_);

    // Reproduction draws on storage too.
    const double fecundity_dt_ = fecundity_dt(mass_storage_, height);
    vars.set_aux(FECUNDITY_DT_ABS_AUX_INDEX, fecundity_dt_);
    vars.set_rate(FECUNDITY_INDEX, fecundity_dt_);
    dmass_storage_dt_ -= fecundity_dt_ * (pars.omega + pars.a_f3);
  } else {
    // Storage is empty: no growth, no reproduction, tissues turn over.
    for (int i : {DHEIGHT_DAREA_LEAF_AUX_INDEX, DAREA_LEAF_DMASS_LIVE_AUX_INDEX,
                  DAREA_LEAF_DT_AUX_INDEX, DMASS_LEAF_DT_AUX_INDEX,
                  DMASS_LEAF_DAREA_LEAF_AUX_INDEX, DMASS_SAPWOOD_DAREA_LEAF_AUX_INDEX,
                  DMASS_SAPWOOD_DT_AUX_INDEX, DMASS_BARK_DAREA_LEAF_AUX_INDEX,
                  DMASS_BARK_DT_AUX_INDEX, DMASS_ROOT_DAREA_LEAF_AUX_INDEX,
                  DMASS_ROOT_DT_AUX_INDEX, FECUNDITY_DT_ABS_AUX_INDEX}) {
      vars.set_aux(i, 0.0);
    }
    vars.set_rate(HEIGHT_INDEX, 0.0);
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_rate(AREA_LEAF_INDEX, -pars.k_l * area_leaf_);
    vars.set_rate(MASS_LEAF_INDEX, -pars.k_l * mass_leaf_);
    vars.set_rate(AREA_HEARTWOOD_INDEX, pars.k_s * area_sapwood_);
    vars.set_rate(MASS_HEARTWOOD_INDEX, pars.k_s * mass_sapwood_);
    vars.set_rate(AREA_SAPWOOD_INDEX, -pars.k_s * area_sapwood_);
    vars.set_rate(MASS_SAPWOOD_INDEX, -pars.k_s * mass_sapwood_);
    vars.set_rate(AREA_BARK_INDEX, -pars.k_b * area_bark_);
    vars.set_rate(MASS_BARK_INDEX, -pars.k_b * mass_bark_);
    vars.set_rate(DIAMETER_STEM_INDEX, 0.0);
    vars.set_rate(MASS_ROOT_INDEX, -pars.k_r * mass_root_);
    vars.set_rate(AREA_STEM_INDEX, -pars.k_b * area_bark_);
  }

  vars.set_rate(MASS_STORAGE_INDEX, dmass_storage_dt_);
  vars.set_rate(MORTALITY_INDEX,
                mortality_dt(storage_portion_, productivity_area_, E_s,
                             vars.state(MORTALITY_INDEX)));
  vars.set_aux(MORTALITY_STORAGE_DEPENDENT_DT_AUX_INDEX,
               mortality_storage_dependent_dt(storage_portion_));
  vars.set_aux(MORTALITY_PRODUCTIVITY_DEPENDENT_DT_AUX_INDEX,
               mortality_productivity_dependent_dt(productivity_area_, E_s));
  vars.set_aux(PRODUCTIVITY_AREA_AUX_INDEX, productivity_area_);
  vars.set_aux(STRESS_AUX_INDEX, E_s);
}

inline void ES20_Model::set_initial_states(const environment_type& environment,
                                           Internals<double>& vars) {
  (void)environment;
  vars.set_state(HEIGHT_INDEX, height_0);
  vars.set_state(AREA_LEAF_INDEX, area_leaf_0);
  vars.set_state(MASS_LEAF_INDEX, mass_leaf_0);
  vars.set_state(AREA_SAPWOOD_INDEX, area_sapwood_0);
  vars.set_state(MASS_SAPWOOD_INDEX, mass_sapwood_0);
  vars.set_state(AREA_BARK_INDEX, area_bark_0);
  vars.set_state(MASS_BARK_INDEX, mass_bark_0);
  vars.set_state(MASS_ROOT_INDEX, mass_root_0);
  vars.set_state(AREA_STEM_INDEX, area_stem_0);
  vars.set_state(MASS_STORAGE_INDEX, mass_storage_0);
  vars.set_state(DIAMETER_STEM_INDEX, diameter_stem_0);
  update_dependent_aux(HEIGHT_INDEX, vars);
  update_dependent_aux(AREA_LEAF_INDEX, vars);
}

// Survival through germination, from the seedling's net production at birth,
// in the environment (and so the stress) of the moment it arrives.
inline double ES20_Model::establishment_probability(
    const environment_type& environment) {
  return establishment_probability_given(net_mass_production_dt(
      environment, height_0, area_leaf_0, mass_leaf_0, mass_sapwood_0,
      mass_bark_0, mass_root_0));
}

inline double ES20_Model::establishment_probability_given(
    double net_mass_production_dt) const {
  if (net_mass_production_dt > 0) {
    const double tmp = pars.a_d0 * area_leaf_0 / net_mass_production_dt;
    return 1.0 / (tmp * tmp + 1.0);
  }
  return 0.0;
}

// The seedling height whose live mass equals the seed mass.
inline double ES20_Model::height_seed() const {
  const double h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
               h1 = height_given_mass_leaf(pars.omega);
  const double tol = this->control.offspring_height_tol;
  auto target = [&](double x) -> double {
    return mass_live_given_height(x) - pars.omega;
  };
  return odelia::root::toms748(target, h0, h1, {tol, tol});
}

inline void ES20_Model::prepare_strategy() {
  function_integrator = quadrature::QK(this->control.function_integration_rule);
  canopy_shape.initialise(pars.eta);
  eta_c = CanopyShape<double>::eta_c(pars.eta);
  crown = canopy_shape.unit_crown(function_integrator);
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
  mass_leaf_0 = mass_leaf(area_leaf_0);
  area_sapwood_0 = area_sapwood(area_leaf_0);
  mass_sapwood_0 = mass_sapwood(area_sapwood_0, height_0);
  area_bark_0 = area_bark(area_leaf_0);
  mass_bark_0 = mass_bark(area_bark_0, height_0);
  mass_root_0 = mass_root(area_leaf_0);
  area_stem_0 = area_stem(area_bark_0, area_sapwood_0, 0.0);
  mass_storage_0 =
      pars.b_s1 * mass_live(mass_leaf_0, mass_bark_0, mass_sapwood_0, mass_root_0);
  diameter_stem_0 = diameter_stem(area_stem_0);
  set_birth_rate_driver();
}

inline ES20_Model::ptr make_strategy_ptr(ES20_Model s) {
  s.prepare_strategy();
  return std::make_shared<ES20_Model>(s);
}

}  // namespace plant

#endif
