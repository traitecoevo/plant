// Built from  inst/include/plant/models/ff16_strategy.h on Wed Sep 16 07:11:11 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_ES20_STRATEGY_H_
#define PLANT_PLANT_ES20_STRATEGY_H_

#include <memory>
#include <plant/strategy.h>
#include <plant/models/es20_environment.h>
#include <plant/models/es20_assimilation.h>

namespace plant {

class ES20_Strategy: public Strategy<ES20_Environment> {
public:
  typedef std::shared_ptr<ES20_Strategy> ptr;
  ES20_Strategy();

  // Overrides ----------------------------------------------

  // update this when the length of state_names changes
  static size_t state_size () { return 15; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
      "area_heartwood",
      "mass_heartwood",
      "mass_storage",
      "area_leaf",
      "mass_leaf",
      "area_sapwood",
      "mass_sapwood",
      "area_bark",
      "mass_bark",
      "mass_root",
      "area_stem",
      "diameter_stem"
      });
  }

  std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "net_mass_production_dt",
      "dbiomass_dt",
      "respiration_dt",
      "respiration_leaf_dt",
      "respiration_sapwood_dt",
      "respiration_root_dt",
      "respiration_bark_dt",
      "darea_leaf_dmass_live",
      "darea_leaf_dt",
      "dmass_leaf_darea_leaf",
      "dmass_sapwood_darea_leaf",
      "dmass_bark_darea_leaf",
      "dmass_root_darea_leaf",
      "dmass_leaf_dt",
      "dmass_sapwood_dt",
      "dmass_bark_dt",
      "dmass_root_dt",
      "dheight_darea_leaf",
      "height_difference",
      "height_adjustment", 
      "mass_sapwood_difference",
      "mass_sapwood_adjustment",
      "fecundity_dt_abs",
      "storage_portion",
      "dead_flag",
      "mortality_storage_dependent_dt",
      "mortality_growth_dependent_dt",
      "productivity_area",
      "stress"
    });
    return ret;
  }

  // Translate generic methods to ES20 strategy leaf area methods

  double competition_effect(double height) const {
    return area_leaf(height);
  }

  /* double competition_effect_state(Internals& vars) const { */
    /* return area_leaf_state(vars); */
  /* } */

  double compute_competition(double z, const Internals& vars) const {
    return area_leaf_above(z, vars.state(HEIGHT_INDEX), vars.state(AREA_LEAF_INDEX));
  }

  void refresh_indices();


  // ES20 Methods  ----------------------------------------------

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

  double mass_above_ground(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_root) const;

  double dbiomass_dt(const ES20_Environment& environment, double mass_storage) const;
  
  double dmass_storage_dt(double net_mass_production_dt_, double dbiomass_dt_ ) const;
  
  double storage_portion(double mass_storage, double mass_leaf, double mass_bark, double mass_sapwood,
                         double mass_root) const;
  
  void compute_rates(const ES20_Environment& environment, bool reuse_intervals,
                Internals& vars);

  void update_dependent_aux(const int index, Internals& vars);


  // [eqn 13] Total maintenance respiration
  double respiration(double mass_leaf, double mass_sapwood,
                     double mass_bark, double mass_root, const ES20_Environment& environment) const;
  
  double respiration_leaf(double mass, const ES20_Environment& environment) const;
  double respiration_bark(double mass, const ES20_Environment& environment) const;
  double respiration_sapwood(double mass, const ES20_Environment& environment) const;
  double respiration_root(double mass, const ES20_Environment& environment) const;

  // [eqn 14] Total turnover
  double turnover(double mass_leaf, double mass_bark,
                  double mass_sapwood, double mass_root) const;
  double turnover_leaf(double mass) const;
  double turnover_bark(double mass) const;
  double turnover_sapwood(double mass) const;
  double turnover_root(double mass) const;

  // [eqn 15] Net production
  double net_mass_production_dt_A(double assimilation, double respiration) const;
  double net_mass_production_dt(const ES20_Environment& environment,
                                double height, double area_leaf, double mass_leaf, double mass_sapwood,
                                double mass_bark, double mass_root,
                                bool reuse_intervals=false);
  
  double net_mass_production_dt(const ES20_Environment& environment,
                                const Internals& vars, 
                                bool reuse_intervals=false);

  // [eqn 17] Rate of offspring production
  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  double darea_leaf_dmass_live(double area_leaf, double height, double mass_sapwood) const;
  
  double mass_leaf_dt(double area_leaf, double darea_leaf_dt_) const; 
  // change in height per change in leaf area
  double dheight_darea_leaf(double area_leaf_, double height_) const;
  // adjustments made to allometry because of turnover time
  double height_area_leaf_allometry_difference(double area_leaf_, double height_) const;
  double height_area_leaf_allometry_adjustment(double area_leaf_, double height_) const;
  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  double dmass_leaf_darea_leaf(double area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  double dmass_sapwood_darea_leaf(double area_leaf, double height, double mass_sapwood) const;
  // adjustments made to allometry because of turnover time
  double mass_sapwood_area_leaf_allometry_difference(double area_leaf_, double height_, double mass_sapwood_) const;
  double mass_sapwood_area_leaf_allometry_adjustment(double area_leaf_, double height_, double mass_sapwood_) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  double dmass_bark_darea_leaf(double area_leaf, double height, double mass_sapwood) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  double dmass_root_darea_leaf() const;
  // Growth rate of basal diameter_stem per unit stem area
  double ddiameter_stem_darea_stem(double area_stem) const;
  // Growth rate of components per unit time:
  double area_leaf_dt(double darea_leaf_dmass_live_, double dbiomass_dt_) const;
  double area_sapwood_dt(double area_leaf_dt) const;
  double mass_sapwood_dt(double dleaf_area_dt_, double dmass_sapwood_darea_leaf_) const;
  double area_heartwood_dt(double area_sapwood) const;
  double area_bark_dt(double area_leaf_dt) const;
  double mass_bark_dt(double dsapwood_mass_dt_) const;
  double area_stem_dt(double area_leaf, double area_leaf_dt) const;
  double diameter_stem_dt(double area_stem, double area_stem_dt) const;
  double mass_root_dt(double area_leaf_dt) const;
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
  
  
  double mortality_dt(double storage_portion, double productivity_area, const ES20_Environment& environment, double cumulative_mortality) const;
  double mortality_growth_independent_dt()const ;
  double mortality_storage_dependent_dt(double storage_portion) const;
  double mortality_growth_dependent_dt(double productivity_area, const ES20_Environment& environment) const;
  // [eqn 20] Survival of seedlings during establishment
  double establishment_probability(const ES20_Environment& environment);
  
  // * Competitive environment
  // [eqn 11] total leaf area above height above height `z` for given plant
  double area_leaf_above(double z, double height, double area_leaf) const;
  // [eqn 10] Fraction of leaf area above height `z`
  double Q(double z, double height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // Set constants within ES20_Strategy
  void prepare_strategy();

  // Initialize all relevant states
  void initialize_states(Internals& vars);

  // * Core traits
  double lma, rho, hmat, omega;
  // * Individual allometry
  // Canopy shape parameters
  double eta, eta_c;
  // Sapwood area per leaf area
  double theta;
  // Empirical constants for scaling relationships
  double a_l1, a_l2, a_r1;
  // Bark area per sapwood area
  double a_b1;
  // * Production
  // Respiration constants
  double r_s, r_b, r_r, r_l;
  // Yield = carbon fixed in tissue per carbon assimilated;
  double a_y;
  // Conversion factor
  double a_bio;
  // Leaf, bark sapwood, and root turnover rates
  double k_l, k_b, k_s, k_r;
  // Leaf productivity parameters  - only used when no N reallocation
  double a_p1, a_p2;

  // * Seed production
  // Accessory cost of reproduction, per seed
  double a_f3;
  // Proportion production allocated to reproduction
  double a_f1;
  // Size range across which individuals mature
  double a_f2;

  // * Mortality
  // Probability of survival during dispersal
  double S_D;
  // Parameter for seedling mortality
  double a_d0;
  // Baseline structural mortality rate
  double d_I;
  // Baseline for storage mortality rate
  double a_dS1;
  // Coefficient for storage production in mortality function
  double a_dS2;
 // Baseline for growth mortality rate
  double a_dG1;
  // Coefficient for dry mass production in mortality function
  double a_dG2;

  // * Storage
  // Time of growth switch (day)
  double t_s;
  // Proportion of storage allocated to biomass
  double a_s;
  // Proportion of storage as related to initial biomass
  double b_s1;
  // Speed of growth switch
  double b_s2;
  
  // Height and leaf area of a (germinated) seed
  double height_0;
  double area_leaf_0;
  double mass_leaf_0;
  double mass_sapwood_0;
  double mass_bark_0;
  double mass_root_0;
  double area_sapwood_0;
  double area_bark_0;
  double area_stem_0;
  double mass_storage_0;
  double diameter_stem_0;

  std::string name;

  ES20_Assimilation assimilator;

  
};

ES20_Strategy::ptr make_strategy_ptr(ES20_Strategy s);

// Built from  src/ff16_strategy.cpp on Wed Sep 16 07:11:11 2020 using the scaffolder, from the strategy:  FF16
#include <iostream>

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
ES20_Strategy::ES20_Strategy() {
  // * Core traits - default values
  lma       = 0.1242302;  // Leaf mass per area [kg / m2]
  rho       = 505.0;      // Wood density [kg/m3]
  hmat      = 30; //16.5958691; // Height at maturation [m]
  omega     = 0.001;//3.8e-5 * 3;     // Seed mass [kg]

  // * Individual allometry
  // Canopy shape parameter (extra calculation here later)
  eta       = 12.0; // [dimensionless]
  // Ratio sapwood area area to leaf area
  theta     = 0.000150887; // [dimensionless]
  // Height - leaf mass scaling
  // a_l1        = 5.44; // height with 1m2 leaf [m]
  // a_l2        = 0.306; // dimensionless scaling of height with leaf area
  // Euc saligna 
  a_l1 = -0.56936;
  a_l2 = 1.48185;
  // Root mass per leaf area
  a_r1        = 0.027691;  //[kg / m]
  // Ratio of bark area : sapwood area
  a_b1         = 0.3768273; // [dimensionless]

  // * Production
  // Ratio of leaf dark respiration to leaf mass [mol CO2 / yr  / kg]
  // =  [mol CO2 / m2 / yr]  |  (39.27 = 2100 * 0.00187)  | narea * photosynthesis_per_nitrogen
  //    / [kg(leaf) / m2 ]   |    / (0.1978791)           | lma
  // Hard coded in value of lma here so that this value doesn't change
  // if that trait changes above.
  r_l   = 39.27 / 0.1978791;
  // Root respiration per mass [mol CO2 / yr / kg]
  r_r   = 217.0;
  // Sapwood respiration per stem mass  [mol CO2 / yr / kg]
  // = respiration per volume [mol CO2 / m3 / yr]
  // /  wood density [kg/m3]
  r_s   = 2.412504;
  // Bark respiration per stem mass
  // assumed to be twice rate of sapwood
  // (NOTE that there is a re-parametrisation here relative to the paper
  // -- r_b is defined (new) as 2*r_s, whereas the paper assumes a
  // fixed multiplication by 2)
  r_b   = 2.0 * r_s;
  // Carbon conversion parameter
  a_y      = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  a_bio  = 2.45e-2;
  // Leaf turnover [/yr]
  k_l=  0.4565855;
 // Bark turnover [/yr]
  k_b    = 0.06;
  // Sapwood turnover [/yr]
  k_s     = 0.06;
  // Root turnover [/yr]
  k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  a_p1   = 151.177775377968; // [mol CO2 / yr / m2]
  a_p2   = 0.204716166503633; // [dimensionless]

  // * Seed production
  // Accessory cost of reproduction
  a_f3  = 3.0 *  3.8e-5; // [kg per seed]

  // Maximum allocation to reproduction
  a_f1   = 0.10; //[dimensionless]
  // Size range across which individuals mature
  a_f2   = 50; // [dimensionless]

// * Mortality parameters
  // Probability of survival during dispersal
  S_D   = 0.25; // [dimensionless]
  // Parameter for seedling survival
  a_d0    = 0.1; //[kg / yr / m2]
  // Baseline for intrinsic mortality
  d_I    = 0.01; // [ / yr]
 // Baseline rate for storage-related mortality
  a_dS1    = 3; // [ / yr]
  // Risk coefficient for dry mass production (per area)
  a_dS2    = 300.0;// [yr m2 / kg ]
  
  // Baseline rate for growth-related mortality
  a_dG1    = 3; // [ / yr]
  // Risk coefficient for dry mass production (per area)
  a_dG2    = 20.0;// [yr m2 / kg ]

  // storage parameters
  // storage allocation parameter [kg kg-1]
  // a_s = 0.06;
  a_s = 0.06 * 365;
  // time of switch [yr]
  t_s = 0.4109589;
  b_s1 = 0.1;
  b_s2 = 120;
  
  // Will get computed properly by prepare_strategy
  height_0 = NA_REAL;
  area_leaf_0 = NA_REAL;
  mass_leaf_0 = NA_REAL;
  mass_sapwood_0 = NA_REAL;
  mass_bark_0 = NA_REAL;
  mass_root_0 = NA_REAL;
  area_sapwood_0 = NA_REAL;
  area_bark_0 = NA_REAL;
  area_stem_0 = NA_REAL;
  mass_storage_0 = NA_REAL;
  eta_c    = NA_REAL;
  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "ES20";
}

void ES20_Strategy::refresh_indices () {
    // Create and fill the name to state index maps
  state_index = std::map<std::string,int>();
  aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (int i = 0; i < state_names_vec.size(); i++) {
    state_index[state_names_vec[i]] = i;
  }
  for (int i = 0; i < aux_names_vec.size(); i++) {
    aux_index[aux_names_vec[i]] = i;
  }
}

// [eqn 2] area_leaf (inverse of [eqn 3])
double ES20_Strategy::area_leaf(double height) const {
  return exp(a_l1 + a_l2 * log(height));
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
double ES20_Strategy::mass_leaf(double area_leaf) const {
  return area_leaf * lma;
}

// [eqn 4] area and mass of sapwood
double ES20_Strategy::area_sapwood(double area_leaf) const {
  return area_leaf * theta;
}

double ES20_Strategy::mass_sapwood(double area_sapwood, double height) const {
  return area_sapwood * height * eta_c * rho;
}

// [eqn 5] area and mass of bark
double ES20_Strategy::area_bark(double area_leaf) const {
  return a_b1 * area_leaf * theta;
}

double ES20_Strategy::mass_bark(double area_bark, double height) const {
  return area_bark * height * eta_c * rho;
}

double ES20_Strategy::area_stem(double area_bark, double area_sapwood,
                            double area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

double ES20_Strategy::diameter_stem(double area_stem) const {
  return std::sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
double ES20_Strategy::mass_root(double area_leaf) const {
  return a_r1 * area_leaf;
}

// [eqn 8] Total mass
double ES20_Strategy::mass_live(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

double ES20_Strategy::mass_total(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_heartwood,
                            double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

double ES20_Strategy::mass_above_ground(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_root;
}

// for updating auxillary state
void ES20_Strategy::update_dependent_aux(const int index, Internals& vars) {
  if (index == AREA_LEAF_INDEX) {
    double area_leaf= vars.state(AREA_LEAF_INDEX);
    vars.set_aux(aux_index.at("competition_effect"), area_leaf);
  }
}


double ES20_Strategy::dbiomass_dt(const ES20_Environment& environment, 
                                  double mass_storage) const {
  if (mass_storage > 0){
    double yr = floor(environment.time);
    double yr_next = yr+1;
    
    // if(time < (yr + (t_s/2))){
    //   return ((a_s)/(1 + exp(-b_s2 * (time - yr))));
    // }
    // else if(time < (yr + t_s + ((1-t_s))/2)){
    //   return ((a_s)/(1 + exp(b_s2 * (time - (t_s+yr)))));
    // }
    // else{
    //   return ((a_s)/(1 + exp(-b_s2 * (time - yr_next))));
    // }
    
    
    if(environment.time < (yr + (t_s/2))){
      return (a_s * mass_storage)/(1 + exp(-b_s2 * (environment.time - yr)));
    }
    else if(environment.time < (yr + t_s + ((1-t_s)/2))){
      return (a_s * mass_storage)/(1 + exp(b_s2 * (environment.time - (t_s+yr))));
    }
    else{
      return (a_s * mass_storage)/(1 + exp(-b_s2 * (environment.time - yr_next)));
    }
  }
  return 0;
}

double ES20_Strategy::dmass_storage_dt(double net_mass_production_dt_, double dbiomass_dt_) const {
  return net_mass_production_dt_ - dbiomass_dt_;
}

double ES20_Strategy::storage_portion(double mass_storage, double mass_leaf, double mass_bark, double mass_sapwood,
                       double mass_root) const {
  return mass_storage / (mass_live(mass_leaf, mass_bark, mass_sapwood, mass_root));
  
}

// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void ES20_Strategy::compute_rates(const ES20_Environment& environment,
                              bool reuse_intervals,
                              Internals& vars) {
  
  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.state(AREA_LEAF_INDEX);
  double area_bark_ = vars.state(AREA_BARK_INDEX);
  double area_heartwood_ = vars.state(AREA_HEARTWOOD_INDEX);
  double area_sapwood_ = vars.state(AREA_SAPWOOD_INDEX);
  double mass_leaf_ = vars.state(MASS_LEAF_INDEX);
  double mass_bark_ = vars.state(MASS_BARK_INDEX);
  double mass_sapwood_ = vars.state(MASS_SAPWOOD_INDEX);
  double mass_root_ = vars.state(MASS_ROOT_INDEX);
  double mass_storage_ = vars.state(MASS_STORAGE_INDEX);
  
  
  
  if(area_leaf_ <= 0){
    area_leaf_ = 0.0001;
    vars.set_state(state_index.at("area_leaf"), area_leaf_);
  }
  
  if(mass_storage_ < 0){
    vars.set_aux(aux_index.at("dead_flag"), 1);
    // mass_storage_ = 0;
    // vars.set_state(state_index.at("mass_storage"), mass_storage_);
  }
  
  // std::cout << "area_leaf: " << area_leaf_ << std::endl;
  // std::cout << "storage 1: " << mass_storage_ << std::endl;
  
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, vars, reuse_intervals);
  
  // std::cout << "A : " << net_mass_production_dt_ << std::endl;
  
  const double dbiomass_dt_ = dbiomass_dt(environment, mass_storage_);
  double dmass_storage_dt_  = dmass_storage_dt(net_mass_production_dt_, dbiomass_dt_);
  const double storage_portion_ = storage_portion(mass_storage_, mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  // 
  // std::cout << "storage 2: " << mass_storage_ << " rate: " << dmass_storage_dt_ << std::endl;
  
  // store the aux sate
  vars.set_aux(aux_index.at("net_mass_production_dt"), net_mass_production_dt_);
  vars.set_aux(aux_index.at("dbiomass_dt"), dbiomass_dt_);
  vars.set_aux(aux_index.at("storage_portion"), storage_portion_);
  vars.set_aux(aux_index.at("respiration_leaf_dt"), a_bio * a_y * respiration_leaf(mass_leaf_,environment));
  vars.set_aux(aux_index.at("respiration_bark_dt"), a_bio * a_y * respiration_bark(mass_bark_,environment));
  vars.set_aux(aux_index.at("respiration_root_dt"), a_bio * a_y * respiration_root(mass_root_,environment));
  vars.set_aux(aux_index.at("respiration_sapwood_dt"), a_bio * a_y * respiration_sapwood(mass_sapwood_,environment));
  vars.set_aux(aux_index.at("respiration_dt"), a_bio * a_y * (respiration_leaf(mass_leaf_,environment) +
    respiration_bark(mass_bark_,environment) + respiration_sapwood(mass_sapwood_,environment) + respiration_root(mass_root_,environment)));
  
  // store the adjustments
  
  vars.set_aux(aux_index.at("height_difference"), height_area_leaf_allometry_difference(area_leaf_, height));
  vars.set_aux(aux_index.at("height_adjustment"), height_area_leaf_allometry_adjustment(area_leaf_, height));
  vars.set_aux(aux_index.at("mass_sapwood_difference"), mass_sapwood_area_leaf_allometry_difference(area_leaf_, height, mass_sapwood_));
  vars.set_aux(aux_index.at("mass_sapwood_adjustment"), mass_sapwood_area_leaf_allometry_adjustment(area_leaf_, height, mass_sapwood_));
  
  // std::cout << "net bio: "<< dbiomass_dt_ << std::endl;
  // 
  if (dbiomass_dt_ > 0) {
    
    // std::cout<< "height " << height << "   leaf area    " << area_leaf_ << std::endl;
    const double dheight_darea_leaf_ = dheight_darea_leaf(area_leaf_, height);
    vars.set_aux(aux_index.at("dheight_darea_leaf"), dheight_darea_leaf_);
    
    // std::cout<< "dheight_darea_leaf " << dheight_darea_leaf_ <<std::endl;
    
    const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_, height, mass_sapwood_);
    const double darea_leaf_dt_ = darea_leaf_dmass_live_ * dbiomass_dt_;
    const double dmass_leaf_dt_ = mass_leaf_dt(area_leaf_, darea_leaf_dt_);
    
    // std::cout<< "ddarea_leaf_dt " << darea_leaf_dt_ <<std::endl;
    
    const double dheight_dt_ = dheight_darea_leaf_ * darea_leaf_dt_;
    
    // std::cout<< "dheight_dt  " << dheight_dt_ <<std::endl;
    
    vars.set_aux(aux_index.at("darea_leaf_dmass_live"), darea_leaf_dmass_live_);
    vars.set_aux(aux_index.at("darea_leaf_dt"), darea_leaf_dt_);
    vars.set_aux(aux_index.at("dmass_leaf_dt"), dmass_leaf_dt_);
    vars.set_aux(aux_index.at("dmass_leaf_darea_leaf"), dmass_leaf_dt_ / darea_leaf_dt_);
    
    vars.set_rate(HEIGHT_INDEX, dheight_dt_);
    vars.set_rate(AREA_LEAF_INDEX, darea_leaf_dt_ - turnover_leaf(area_leaf_));
    vars.set_rate(MASS_LEAF_INDEX, dmass_leaf_dt_ - turnover_leaf(mass_leaf_));
    
    // Changes in sapwood and heartwood
    
    const double darea_heartwood_dt_ = area_heartwood_dt(area_sapwood_);
    const double dmass_heartwood_dt_ = mass_heartwood_dt(mass_sapwood_);
    
    vars.set_rate(state_index.at("area_heartwood"), darea_heartwood_dt_);
    vars.set_rate(state_index.at("mass_heartwood"), dmass_heartwood_dt_);
    
    const double darea_sapwood_dt_ = area_sapwood_dt(darea_leaf_dt_);
    const double dmass_sapwood_darea_leaf_ = dmass_sapwood_darea_leaf(area_leaf_, height, mass_sapwood_);
    const double dmass_sapwood_dt_ = mass_sapwood_dt(darea_leaf_dt_, dmass_sapwood_darea_leaf_);
    
    vars.set_aux(aux_index.at("dmass_sapwood_darea_leaf"), dmass_sapwood_darea_leaf_);
    vars.set_aux(aux_index.at("dmass_sapwood_dt"), dmass_sapwood_dt_);
    
    vars.set_rate(state_index.at("area_sapwood"), darea_sapwood_dt_ - turnover_sapwood(area_sapwood_));
    vars.set_rate(state_index.at("mass_sapwood"), dmass_sapwood_dt_ - turnover_sapwood(mass_sapwood_));
    
    // changes in bark
    const double darea_bark_dt_ = area_bark_dt(darea_leaf_dt_);
    const double dmass_bark_dt_ = mass_bark_dt(dmass_sapwood_dt_);
    
    vars.set_aux(aux_index.at("dmass_bark_darea_leaf"), dmass_bark_dt_ / darea_leaf_dt_);
    vars.set_aux(aux_index.at("dmass_bark_dt"), dmass_bark_dt_);
    
    vars.set_rate(state_index.at("area_bark"), darea_bark_dt_ - turnover_bark(area_bark_));
    vars.set_rate(state_index.at("mass_bark"), dmass_bark_dt_ - turnover_bark(mass_bark_));
    
    // changes in stem
    const double darea_stem_dt_ = darea_sapwood_dt_ + darea_heartwood_dt_ + darea_bark_dt_ - turnover_bark(area_bark_);
    const double ddiameter_stem_dt_ = diameter_stem_dt(area_stem(area_heartwood_, area_sapwood_, area_bark_), darea_stem_dt_);
    
    vars.set_rate(state_index.at("diameter_stem"), ddiameter_stem_dt_);
    vars.set_rate(state_index.at("area_stem"), darea_stem_dt_);
    
    // changes in root
    const double dmass_root_dt_ = mass_root_dt(darea_leaf_dt_);
    
    vars.set_aux(aux_index.at("dmass_root_darea_leaf"), dmass_root_dt_ / darea_leaf_dt_);
    vars.set_aux(aux_index.at("dmass_root_dt"), dmass_root_dt_);
    
    vars.set_rate(state_index.at("mass_root"), dmass_root_dt_ - turnover_root(mass_root_));
    
    
    vars.set_aux(aux_index.at("fecundity_dt_abs"), fecundity_dt(mass_storage_, height));
    
    // fecundity 
    vars.set_rate(FECUNDITY_INDEX,
                  fecundity_dt(mass_storage_, height));
    
    // update losses due to reproduction
    dmass_storage_dt_ = dmass_storage_dt_ - fecundity_dt(mass_storage_, height) * (omega + a_f3);
    
  } else {
    vars.set_rate(HEIGHT_INDEX, 0.0);
    // std::cout<< "Height: " << height << " Rate: " << 0 << std::endl;
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_aux(aux_index.at("fecundity_dt_abs"), 0);
    vars.set_rate(AREA_LEAF_INDEX, - turnover_leaf(area_leaf_));
    vars.set_rate(MASS_LEAF_INDEX, - turnover_leaf(mass_leaf_));
    vars.set_rate(state_index.at("area_heartwood"), turnover_sapwood(area_sapwood_));
    vars.set_rate(state_index.at("mass_heartwood"), turnover_sapwood(mass_sapwood_));
    vars.set_rate(state_index.at("area_sapwood"), - turnover_sapwood(area_sapwood_));
    vars.set_rate(state_index.at("mass_sapwood"), (- turnover_sapwood(mass_sapwood_))); // <- error in this (?)
    vars.set_rate(state_index.at("area_bark"), - turnover_bark(area_bark_));
    vars.set_rate(state_index.at("mass_bark"), - turnover_bark(mass_bark_));
    vars.set_rate(state_index.at("diameter_stem"), 0.0);
    vars.set_rate(state_index.at("mass_root"), - turnover_root(mass_root_));
    vars.set_rate(state_index.at("area_stem"), - turnover_bark(area_bark_));
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MASS_STORAGE_INDEX, dmass_storage_dt_);
  vars.set_rate(MORTALITY_INDEX,
                mortality_dt(storage_portion_, net_mass_production_dt_ / area_leaf_, environment, vars.state(MORTALITY_INDEX)));
  
  vars.set_aux(aux_index.at("mortality_storage_dependent_dt"), mortality_storage_dependent_dt(storage_portion_));
  vars.set_aux(aux_index.at("mortality_growth_dependent_dt"), mortality_growth_dependent_dt(net_mass_production_dt_ / area_leaf_, environment));
  vars.set_aux(aux_index.at("productivity_area"), net_mass_production_dt_ / area_leaf_);
  vars.set_aux(aux_index.at("stress"), environment.getStress());
  
  
  // std::cout<< "Finished compute rates "<< std::endl;
}


// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by a_y*a_bio.
// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by a_y*a_bio.
double ES20_Strategy::respiration(double mass_leaf, double mass_sapwood,
                                  double mass_bark, double mass_root, const ES20_Environment& environment) const {
  return respiration_leaf(mass_leaf, environment) +
    respiration_bark(mass_bark,environment) +
    respiration_sapwood(mass_sapwood,environment) +
    respiration_root(mass_root, environment);
}

double ES20_Strategy::respiration_leaf(double mass, const ES20_Environment& environment) const {
  const double stress = environment.getStress();
  // std::cout<< "Stress respiration : " << (0.5 + stress/2) ;
  return r_l * mass * (0.5 + stress/2);
}

double ES20_Strategy::respiration_bark(double mass, const ES20_Environment& environment) const {
  const double stress = environment.getStress();
  return r_b * mass * (0.5 + stress/2);
}

double ES20_Strategy::respiration_sapwood(double mass, const ES20_Environment& environment) const {
  const double stress = environment.getStress();
  return r_s * mass * (0.5 + stress/2);
}

double ES20_Strategy::respiration_root(double mass, const ES20_Environment& environment) const {
  const double stress = environment.getStress();
  return r_r * mass * (0.5 + stress/2);
}


// [eqn 14] Total turnover
double ES20_Strategy::turnover(double mass_leaf, double mass_bark,
                          double mass_sapwood, double mass_root) const {
   return turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root);
}

double ES20_Strategy::turnover_leaf(double mass) const {
  return k_l * mass;
}

double ES20_Strategy::turnover_bark(double mass) const {
  return k_b * mass;
}

double ES20_Strategy::turnover_sapwood(double mass) const {
  return k_s * mass;
}

double ES20_Strategy::turnover_root(double mass) const {
  return k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
double ES20_Strategy::net_mass_production_dt_A(double assimilation, double respiration) const {
  return a_bio * a_y * (assimilation - respiration);
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
double ES20_Strategy::net_mass_production_dt(const ES20_Environment& environment,
                                             double height, double area_leaf, double mass_leaf, double mass_sapwood,
                                             double mass_bark, double mass_root,
                                             bool reuse_intervals) {
  const double assimilation_ = assimilator.assimilate(control, environment, height,
                                                      area_leaf, reuse_intervals);
  const double respiration_ =
    respiration(mass_leaf, mass_sapwood, mass_bark, mass_root, environment);
  return net_mass_production_dt_A(assimilation_, respiration_);
}

double ES20_Strategy::net_mass_production_dt(const ES20_Environment& environment, 
                                             const Internals& vars, 
                                             bool reuse_intervals) {
  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.state(AREA_LEAF_INDEX);
  double mass_leaf_ = vars.state(MASS_LEAF_INDEX);
  double mass_bark_ = vars.state(MASS_BARK_INDEX);
  double mass_sapwood_ = vars.state(MASS_SAPWOOD_INDEX);
  double mass_root_ = vars.state(MASS_ROOT_INDEX);
  
  return(net_mass_production_dt(environment, height, area_leaf_, mass_leaf_, mass_sapwood_, mass_bark_, mass_root_, reuse_intervals));
}

// [eqn 17] Rate of offspring production
double ES20_Strategy::fecundity_dt(double net_mass_production_dt,
                                   double fraction_allocation_reproduction) const {
  double mass_storage = net_mass_production_dt;
  double height = fraction_allocation_reproduction;
  return mass_storage * a_f1 / ((omega + a_f3) * (1.0 + exp(a_f2 * (1.0 - height / hmat))));
}

double ES20_Strategy::area_leaf_dt(double darea_leaf_dmass_live_, double dbiomass_dt_) const{
  return darea_leaf_dmass_live_ * dbiomass_dt_;
}

double ES20_Strategy::mass_leaf_dt(double area_leaf, double darea_leaf_dt_) const {
  return dmass_leaf_darea_leaf(area_leaf) * darea_leaf_dt_ ;
}

double ES20_Strategy::darea_leaf_dmass_live(double area_leaf, double height, double mass_sapwood) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
                  + dmass_sapwood_darea_leaf(area_leaf, height, mass_sapwood)
                  + dmass_bark_darea_leaf(area_leaf, height, mass_sapwood)
                  + dmass_root_darea_leaf());
}

// TODO: Ordering below here needs working on, probably as @dfalster
// does equation documentation?
double ES20_Strategy::dheight_darea_leaf(double area_leaf_, double height_) const {
  const double adjustment = height_area_leaf_allometry_adjustment(area_leaf_, height_);
  return (exp(log(area_leaf_)/a_l2 - a_l1/a_l2) /(area_leaf_ * a_l2)) * adjustment;
}

double ES20_Strategy::height_area_leaf_allometry_difference(double area_leaf_, double height_) const {
  return (area_leaf_ - area_leaf(height_));
}

double ES20_Strategy::height_area_leaf_allometry_adjustment(double area_leaf_, double height_) const {
  return 2 / (1 + exp(-40 * (height_area_leaf_allometry_difference(area_leaf_, height_))));
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
double ES20_Strategy::dmass_leaf_darea_leaf(double /* area_leaf */) const {
  return lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
double ES20_Strategy::dmass_sapwood_darea_leaf(double area_leaf_, double height_, double mass_sapwood_) const {
  const double adjustment = mass_sapwood_area_leaf_allometry_adjustment(area_leaf_, height_, mass_sapwood_);
  return rho * eta_c *theta * (height_ + area_leaf_ * dheight_darea_leaf(area_leaf_, height_)) * adjustment;
}

double ES20_Strategy::mass_sapwood_area_leaf_allometry_difference(double area_leaf_, double height_, double mass_sapwood_) const {
  const double difference = mass_sapwood(area_sapwood(area_leaf_), height_) - mass_sapwood_;
  return difference;
}

double ES20_Strategy::mass_sapwood_area_leaf_allometry_adjustment(double area_leaf_, double height_, double mass_sapwood_) const {
  const double adjustment = 2 / (1 + exp(-40 * (mass_sapwood_area_leaf_allometry_difference(area_leaf_, height_, mass_sapwood_))));
  return adjustment;
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
double ES20_Strategy::dmass_bark_darea_leaf(double area_leaf, double height, double mass_sapwood) const {
  return a_b1 * dmass_sapwood_darea_leaf(area_leaf, height, mass_sapwood);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
double ES20_Strategy::dmass_root_darea_leaf() const {
  return a_r1;
}

// Growth rate of basal diameter_stem per unit time
double ES20_Strategy::ddiameter_stem_darea_stem(double area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
double ES20_Strategy::area_sapwood_dt(double area_leaf_dt) const {
  return area_leaf_dt * theta;
}

double ES20_Strategy::mass_sapwood_dt(double dleaf_area_dt_, double dmass_sapwood_darea_leaf_) const{
  return dleaf_area_dt_ * dmass_sapwood_darea_leaf_;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
double ES20_Strategy::area_heartwood_dt(double area_leaf) const {
  return k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
double ES20_Strategy::area_bark_dt(double area_leaf_dt) const {
  return a_b1 * area_leaf_dt * theta;
}

// Growth rate of bark mass at base per unit time
double ES20_Strategy::mass_bark_dt(double dsapwood_mass_dt) const {
  return a_b1 * dsapwood_mass_dt;
}

// Growth rate of stem basal area per unit time
double ES20_Strategy::area_stem_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
double ES20_Strategy::diameter_stem_dt(double area_stem, double area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
double ES20_Strategy::mass_root_dt(double area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf();
}

double ES20_Strategy::mass_live_dt(double fraction_allocation_reproduction,
                               double net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

// TODO: Change top two to use mass_live_dt
double ES20_Strategy::mass_total_dt(double fraction_allocation_reproduction,
                                     double net_mass_production_dt,
                                     double mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
double ES20_Strategy::mass_above_ground_dt(double area_leaf,
                                       double fraction_allocation_reproduction,
                                       double net_mass_production_dt,
                                       double mass_heartwood_dt,
                                       double area_leaf_dt) const {
  const double mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf();
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

double ES20_Strategy::mass_heartwood_dt(double mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


double ES20_Strategy::mass_live_given_height(double height) const {
  double area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

double ES20_Strategy::height_given_mass_leaf(double mass_leaf) const {
  return exp(1/a_l2 * log(mass_leaf / lma) - a_l1/a_l2);
}

double ES20_Strategy::mortality_dt(double storage_portion,
                                   double productivity_area,
                                   const ES20_Environment& environment,
                                   double cumulative_mortality) const {
  
  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (R_FINITE(cumulative_mortality)) {
    return
      mortality_growth_independent_dt() +
      mortality_storage_dependent_dt(storage_portion) +
      mortality_growth_dependent_dt(productivity_area, environment);
  } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    return 0.0;
  }
}

double ES20_Strategy::mortality_growth_independent_dt() const {
  return d_I;
}

double ES20_Strategy::mortality_storage_dependent_dt(double storage_portion) const {
  return a_dS1 * exp(-a_dS2 * storage_portion);
}

double ES20_Strategy::mortality_growth_dependent_dt(double productivity_area, const ES20_Environment& environment) const{
  //when stressed mortality is purely dependent on storage 
  const double stress = environment.getStress();
  return stress * a_dG1 * exp(-a_dG2 * productivity_area);
}

// [eqn 20] Survival of seedlings during establishment
double ES20_Strategy::establishment_probability(const ES20_Environment& environment) {
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0, mass_leaf_0, mass_sapwood_0,
                           mass_bark_0, mass_root_0);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0);
  } else {
    return 0.0;
  }
}

double ES20_Strategy::area_leaf_above(double z, double height, double area_leaf) const {
  // std::cout << "In Strategy . " << "Leaf area " << area_leaf << "   Q    " <<  Q(z, height) << "   final values:   " << area_leaf * Q(z, height) <<std::endl;
  return area_leaf * Q(z, height);
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height'
double ES20_Strategy::Q(double z, double height) const {
  if (z > height) {
    return 0.0;
  }
  const double tmp = 1.0-pow(z / height, eta);
  return tmp * tmp;
}

// The aim is to find a plant height that gives the correct seed mass.
double ES20_Strategy::height_seed(void) const {
  
  
  // If height is preset then use that as a start
  if(height_0 == height_0){
    return height_0;
  }

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const double
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(omega);

  const double tol = control.plant_seed_tol;
  const size_t max_iterations = control.plant_seed_iterations;

  auto target = [&] (double x) mutable -> double {
    return mass_live_given_height(x) - omega;
  };

  return util::uniroot(target, h0, h1, tol, max_iterations);
}

void ES20_Strategy::prepare_strategy() {
  // Set up the integrator
  control.initialize();
  assimilator.initialize(a_p1, a_p2, eta);
  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2/(1 + eta) + 1/(1 + 2*eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
  mass_leaf_0 = mass_leaf(area_leaf_0);
  area_sapwood_0 = area_sapwood(area_leaf_0);
  mass_sapwood_0 = mass_sapwood(area_sapwood_0, height_0);
  area_bark_0 = area_bark(area_leaf_0);
  mass_bark_0 = mass_bark(area_bark_0, height_0);
  mass_root_0 = mass_root(area_leaf_0);
  area_stem_0 = area_stem(area_bark_0, area_sapwood_0, 0); 
  mass_storage_0 = b_s1 * mass_live(mass_leaf_0, mass_bark_0, mass_sapwood_0, mass_root_0);
  diameter_stem_0 = diameter_stem(area_stem_0);
}

void ES20_Strategy::initialize_states(Internals &vars){
  vars.set_state(HEIGHT_INDEX, height_0);
  vars.set_state(AREA_LEAF_INDEX, area_leaf_0);
  vars.set_state(MASS_LEAF_INDEX, mass_leaf_0);
  vars.set_state(state_index.at("area_sapwood"), area_sapwood_0);
  vars.set_state(state_index.at("mass_sapwood"), mass_sapwood_0);
  vars.set_state(state_index.at("area_bark"), area_bark_0);
  vars.set_state(state_index.at("mass_bark"), mass_bark_0);
  vars.set_state(state_index.at("mass_root"), mass_root_0);
  vars.set_state(state_index.at("area_stem"), area_stem_0);
  vars.set_state(state_index.at("mass_storage"), mass_storage_0);
  vars.set_state(state_index.at("diameter_stem"), diameter_stem_0);
  update_dependent_aux(AREA_LEAF_INDEX, vars);
  
}

ES20_Strategy::ptr make_strategy_ptr(ES20_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<ES20_Strategy>(s);
}

}

#endif
