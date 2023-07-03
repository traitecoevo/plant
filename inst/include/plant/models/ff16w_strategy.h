// Built from  inst/include/plant/models/ff16_strategy.h using the scaffolder,
// from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_WATER_STRATEGY_H_
#define PLANT_PLANT_WATER_STRATEGY_H_

#include <plant/leaf_model.h>
#include <plant/models/ff16_strategy.h>

namespace plant {

class FF16w_Strategy : public FF16_Strategy {
public:
  typedef std::shared_ptr<FF16w_Strategy> ptr;
  FF16w_Strategy();


  double compute_average_light_environment(double z, double height,
                                           const FF16_Environment &environment);

  // calculate the amount of water transpired relativised by leaf area index.

  double evapotranspiration_dt(double area_leaf_);

  // Net production

 double darea_leaf_dmass_live(double area_leaf) const;

virtual double net_mass_production_dt(const FF16_Environment &environment,
                                        double height, double area_leaf_, Internals &vars, 
                                        bool reuse_intervals = false);
                                        
  virtual void compute_rates(const FF16_Environment &environment,
                             bool reuse_intervals, Internals &vars);
                             
  virtual void prepare_strategy();

// names of auxillary variables
virtual std::vector<std::string> aux_names() {
  std::vector<std::string> ret(
      {"competition_effect", "net_mass_production_dt", "opt_psi_stem_","opt_ci_", "count","profit_","assim_colimited_","hydraulic_cost_", "ci_","darea_leaf_dmass_live_", "stom_cond_CO2_", "transpiration_"});
  // add the associated computation to compute_rates and compute there
  if (collect_all_auxiliary) {
    ret.push_back("area_sapwood");
  }
  return ret;
}

  Leaf leaf;

  // leaf traits - default values Eucalyptus saligna
  double vcmax_25 = 96;
  double p_50 = 3.4;
  double K_s = (pow(p_50/1.731347,1/-0.7246377))*2;
  double c = log(log(1-0.5)/log(1-0.88))/(log(p_50) - log(5.16));
  double b = p_50 / std::pow(-log(1 - 50.0 / 100.0), 1 / c);
  double psi_crit = b*std::pow(log(1/0.05),1/c); // derived from b and c
  double epsilon_leaf = 0.001;
  double beta1 = 20000;
  double beta2 = 1.5;
  double jmax_25 = vcmax_25*1.64;
  double hk_s = 4;
  double a = 0.30; // effective quantum yield of electron transport  (mol photon mol ^-1 electron)  Sabot et al. 2020
  double curv_fact_elec_trans = 0.7; 
  double curv_fact_colim = 0.99; 
  double B_rs1 = 0.01;
  double B_lf2 = 0.01;
  double B_lf3 = 0.01;  
  double B_lf5 = 0.01;

  //nitrogen allocation traits (parameterised from Austraits 4.1.0)
  double nmass_l = 13e-3; // kg N kg^-1 mass
  double nmass_s = 1.98e-3; // kg N kg^-1 mass
  double nmass_b = 3.40e-3; // kg N kg^-1 mass
  double nmass_r = 3.35e-3; // kg N kg^-1 mass
  double dmass_dN = 0; //change in mass per change in kg kg^-1 N
};


FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s);

} // namespace plant

#endif
