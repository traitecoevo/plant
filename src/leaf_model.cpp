#include <plant/leaf_model.h>
#include <cmath>

namespace plant {
Leaf::Leaf()
    :
    vcmax_25(96), // umol m^-2 s^-1 
    c(2.680147), //unitless
    b(3.898245), //-MPa
    psi_crit(5.870283), //-MPa 
    beta1(20000), //hydraulic cost for Bartlett method umol m^-3 s^-1
    beta2(1.5), //exponent for effect of hydraulic risk (unitless)
    jmax_25(157.44), // maximum electron transport rate umol m^-2 s^-1
    hk_s(4),  // maximum hydraulic-dependent sapwood turnover rate yr ^ -1
    a(0.30), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(0.7), //curvature factor for the light response curve (unitless)
    curv_fact_colim(0.99), //curvature factor for the colimited photosythnthesis equatiom
    GSS_tol_abs(1e-3),
    vulnerability_curve_ncontrol(100),
    ci_abs_tol(1e-3),
    ci_niter(1000)
   {
      setup_transpiration(100); // arg: num control points for integration
      setup_clean_leaf();
}

Leaf::Leaf(double vcmax_25, double c, double b,
           double psi_crit, // derived from b and c,
           double beta1, double beta2, double jmax_25, double hk_s,
           double a, double curv_fact_elec_trans, double curv_fact_colim, 
           double GSS_tol_abs,
           double vulnerability_curve_ncontrol,
           double ci_abs_tol,
           double ci_niter)
    : vcmax_25(vcmax_25), // umol m^-2 s^-1 
    c(c), //unitless
    b(b), //-MPa
    psi_crit(psi_crit), //-MPa 
    beta1(beta1), //hydraulic cost for Bartlett method umol m^-3 s^-1
    beta2(beta2), //exponent for effect of hydraulic risk (unitless)
    jmax_25(jmax_25), // maximum electron transport rate umol m^-2 s^-1
    hk_s(hk_s),  // maximum hydraulic-dependent sapwood turnover rate yr ^ -1
    a(a), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(curv_fact_elec_trans), //curvature factor for the light response curve (unitless)
    curv_fact_colim(curv_fact_colim), //curvature factor for the colimited photosythnthesis equation
    GSS_tol_abs(GSS_tol_abs),
    vulnerability_curve_ncontrol(vulnerability_curve_ncontrol),
    ci_abs_tol(ci_abs_tol),
    ci_niter(ci_niter)
   {
      setup_transpiration(vulnerability_curve_ncontrol); // arg: num control points for integration
      setup_clean_leaf();
}


//set various states and physiology parameters obtained from ff16w to NA to clean leaf object
void Leaf::setup_clean_leaf() {
  ci_ = NA_REAL; // Pa
  stom_cond_CO2_= NA_REAL; //mol Co2 m^-2 s^-1 
  assim_colimited_= NA_REAL; // umol C m^-2 s^-1 
  transpiration_= NA_REAL; // kg m^-2 s^-1 
  profit_= NA_REAL; // umol C m^-2 s^-1 
  lambda_= NA_REAL; // umol C m^-2 s^-1 kg^-1 m^2 s^1
  lambda_analytical_= NA_REAL; // umol C m^-2 s^-1 kg^-1 m^2 s^1
  hydraulic_cost_= NA_REAL; // umol C m^-2 s^-1 
  electron_transport_= NA_REAL; //electron transport rate umol m^-2 s^-1
  gamma_= NA_REAL;
  ko_= NA_REAL;
  kc_= NA_REAL;
  km_= NA_REAL;
  R_d_= NA_REAL;
  leaf_specific_conductance_max_= NA_REAL; //kg m^-2 s^-1 MPa^-1 
  sapwood_volume_per_leaf_area_ = NA_REAL; //m^3 SA m^-2 LA
  rho_= NA_REAL; //kg m^-3
  vcmax_= NA_REAL; //kg m^-3
  jmax_= NA_REAL; //kg m^-3
  a_bio_= NA_REAL; //kg mol^-1
  psi_soil_= NA_REAL; //-MPa
  leaf_temp_= NA_REAL; // deg C
  PPFD_= NA_REAL; //umol m^-2 s^-1
  atm_vpd_= NA_REAL; //kPa 
  atm_o2_kpa_= NA_REAL; // kPa
  atm_kpa_= NA_REAL; // kPa
  ca_= NA_REAL; //Pa
  opt_psi_stem_= NA_REAL; //-MPa 
  opt_ci_= NA_REAL; //Pa 
}

//sets various parameters which are constant for a given node at a given time

void Leaf::set_physiology(double rho, double a_bio, double PPFD, double psi_soil, double leaf_specific_conductance_max, double atm_vpd, double ca, double sapwood_volume_per_leaf_area, double leaf_temp, double atm_o2_kpa, double atm_kpa) {
   rho_ = rho;
   a_bio_ = a_bio;
   atm_vpd_ = atm_vpd;
   leaf_temp_ = leaf_temp;
   atm_kpa_ = atm_kpa;
   atm_o2_kpa_ = atm_o2_kpa;
   PPFD_ = PPFD;
   psi_soil_ = psi_soil;
   leaf_specific_conductance_max_ = leaf_specific_conductance_max;
   sapwood_volume_per_leaf_area_ = sapwood_volume_per_leaf_area;
   ca_ = ca;
   vcmax_ = peak_arrh_curve(vcmax_ha, vcmax_25, leaf_temp_, vcmax_H_d, vcmax_d_S);
   jmax_ = peak_arrh_curve(jmax_ha, jmax_25, leaf_temp_, jmax_H_d, jmax_d_S);
   electron_transport_ = electron_transport();
   gamma_ = arrh_curve(gamma_ha, gamma_25, leaf_temp_);
   ko_ = arrh_curve(ko_ha, ko_25, leaf_temp_);
   kc_ = arrh_curve(kc_ha, kc_25, leaf_temp_);
   R_d_ = vcmax_*0.015;
   km_ = (kc_*umol_per_mol_to_Pa)*(1 + (atm_o2_kpa_*kPa_to_Pa)/(ko_*umol_per_mol_to_Pa));


// set lambda, if psi_soil is higher than psi_crit, then set to 0. Currently doing both the numerical and analytical version. Ideally would do one.
  // if(psi_soil >= psi_crit){
  //   lambda_ = 0;
  //   lambda_analytical_ = 0;

  // } else {
  //   set_leaf_states_rates_from_psi_stem(psi_crit);
  //   lambda_ = assim_colimited(ci_) / hydraulic_cost_Sperry(psi_crit);
  //   set_leaf_states_rates_from_psi_stem_analytical(psi_crit);
  //   lambda_analytical_ = assim_colimited_analytical(ci_) / hydraulic_cost_Sperry(psi_crit);
  // }
}

double Leaf::arrh_curve(double Ea, double ref_value, double leaf_temp) const {
  return ref_value*exp(Ea*((leaf_temp+C_to_K) - (25 + C_to_K))/((25 + C_to_K)*R*(leaf_temp+C_to_K)));
}

double Leaf::peak_arrh_curve(double Ea, double ref_value, double leaf_temp, double H_d, double d_S) const {
  double arrh = arrh_curve(Ea, ref_value, leaf_temp);
  double arg2 = 1 + exp((d_S*(25 + C_to_K) - H_d)/(R*(25 + C_to_K)));
  double arg3 = 1 + exp((d_S*(leaf_temp + C_to_K) - H_d)/(R*(leaf_temp + C_to_K)));

  return arrh * arg2/arg3;
}


// transpiration supply functions

// returns proportion of conductance taken from hydraulic vulnerability curve (unitless)
double Leaf::proportion_of_conductivity(double psi) const {
  return exp(-pow((psi / b), c));
}

// set spline for proportion of conductivity
void Leaf::setup_transpiration(double resolution) {
  // integrate and accumulate results
  auto x_psi_ = std::vector<double>{0.0};  // {0.0}
  auto y_cumulative_transpiration_ = std::vector<double>{0.0}; // {0.0}
  double step = (b*pow((log(1/0.01)),(1/c)))/resolution;
  
  for (double psi_spline = 0.0 + step; psi_spline <= (b*pow((log(1/0.01)),(1/c))); psi_spline += step) {

    double E_psi = step * ((proportion_of_conductivity(psi_spline-step) + proportion_of_conductivity(psi_spline))/2) + y_cumulative_transpiration_.back();
    x_psi_.push_back(psi_spline); // x values for spline
    y_cumulative_transpiration_.push_back(E_psi); // y values for spline
}
// setup interpolator
transpiration_from_psi.init(x_psi_, y_cumulative_transpiration_);
transpiration_from_psi.set_extrapolate(false);

psi_from_transpiration.init(y_cumulative_transpiration_, x_psi_);
psi_from_transpiration.set_extrapolate(false);
}

// replace f with some other function, returns E kg m^-2 s^-1

double Leaf::transpiration_full_integration(double psi_stem) {
  std::function<double(double)> f;
  f = [&](double psi) -> double { return proportion_of_conductivity(psi); };
  
  return leaf_specific_conductance_max_ * integrator.integrate(f, psi_soil_, psi_stem);
 }

//calculates supply-side transpiration from psi_stem and psi_soil, returns kg h20 s^-1 m^-2 LA
double Leaf::transpiration(double psi_stem) {
  // integration of proportion_of_conductivity over [psi_soil_, psi_stem]
  return leaf_specific_conductance_max_ * (transpiration_from_psi.eval(psi_stem) - transpiration_from_psi.eval(psi_soil_));
  // return (transpiration_full_integration(psi_stem));

  
}

// converts a known transpiration to its corresponding psi_stem, returns -MPa
double Leaf::transpiration_to_psi_stem(double transpiration_) {
  // integration of proportion_of_conductivity over [psi_soil_, psi_stem]
  double E_psi_stem = transpiration_/leaf_specific_conductance_max_ +  transpiration_from_psi.eval(psi_soil_);
  return psi_from_transpiration.eval(E_psi_stem);
  }

// returns stomatal conductance to CO2, mol C m^-2 LA s^-1
double Leaf:: stom_cond_CO2(double psi_stem) {
  double transpiration_ = transpiration(psi_stem);
  return atm_kpa_ * transpiration_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
}


// biochemical photosynthesis model equations
//ensure that units of PPFD_ actually correspond to something real.
// electron trnansport rate based on light availability and vcmax assuming co-limitation hypothesis
double Leaf::electron_transport() {



  double electron_transport_ = (a * PPFD_ + jmax_ - sqrt(pow(a * PPFD_ + jmax_, 2) - 
  4 * curv_fact_elec_trans * a * PPFD_ * jmax_)) / (2 * curv_fact_elec_trans); // check brackets are correct

  // double electron_transport_ = (4*a*PPFD_)/sqrt(pow(4*a*PPFD_/jmax_,2)+ 1);
    return electron_transport_;           
}

//calculate the rubisco-limited assimilation rate, returns umol m^-2 s^-1
double Leaf::assim_rubisco_limited(double ci_) {

  return (vcmax_ * (ci_ - gamma_ * umol_per_mol_to_Pa)) / (ci_ + km_);

}

//calculate the light-limited assimilation rate, returns umol m^-2 s^-1
double Leaf::assim_electron_limited(double ci_) {
  

  return electron_transport_ / 4 *
  ((ci_ - gamma_ * umol_per_mol_to_Pa) / (ci_ + 2 * gamma_ * umol_per_mol_to_Pa));
}

// returns co-limited assimilation umol m^-2 s^-1
double Leaf::assim_colimited(double ci_) {
  
  double assim_rubisco_limited_ = assim_rubisco_limited(ci_) ;
  double assim_electron_limited_ = assim_electron_limited(ci_);

  // no dark respiration included at the moment
  return (assim_rubisco_limited_ + assim_electron_limited_ - sqrt(pow(assim_rubisco_limited_ + assim_electron_limited_, 2) - 4 * curv_fact_colim * assim_rubisco_limited_ * assim_electron_limited_)) /
             (2 * curv_fact_colim)- R_d_;


}

// returns co-limited assimilation based only on light-limited transport rate and empirically-parameterised smothing term, returns umol m^-2 s^-1
double Leaf::assim_colimited_analytical(double ci_) {

  double c2 = 13.13652;

  // no dark respiration included at the moment
  
  return electron_transport_ / 4 *
         ((ci_ - gamma_ * umol_per_mol_to_Pa) /
          (ci_ + c2));
}

// A - gc curves

// returns difference between co-limited assimilation and stom_cond_CO2, to be minimised (umol m^-2 s^-1)
double Leaf::assim_minus_stom_cond_CO2(double x, double psi_stem) {

  double assim_colimited_x_ = assim_colimited(x) + R_d_;

  double stom_cond_CO2_x_ = stom_cond_CO2(psi_stem);
  return assim_colimited_x_ * umol_to_mol -
         (stom_cond_CO2_x_ * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
}

// converts psi stem to ci, used to find ci which makes A(ci) = gc(ca - ci)
double Leaf::psi_stem_to_ci(double psi_stem) {
  // not clear what x is here
  

  auto target = [&](double x) mutable -> double {
    return assim_minus_stom_cond_CO2(x, psi_stem);
  };

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  return ci_ = util::uniroot(target, gamma_ * umol_per_mol_to_Pa, ca_, ci_abs_tol, ci_niter);
}

// given psi_stem, find assimilation, transpiration and stomal conductance to c02
void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem) {
  
  if (psi_soil_ >= psi_stem){
    ci_ = gamma_*umol_per_mol_to_Pa;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    } else{
      ci_ = psi_stem_to_ci(psi_stem);
      transpiration_ = transpiration(psi_stem);
      stom_cond_CO2_ = atm_kpa_ * transpiration_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
      }
  
  assim_colimited_ = assim_colimited(ci_);
  

}

//this set of equations can make values greater than 40 (ca) when PPFD is extremely low
double Leaf::psi_stem_to_ci_analytical(double psi_stem) {

    stom_cond_CO2_ = stom_cond_CO2(psi_stem);
    double c2 = 13.13652;
        
    double first_term = 8 * electron_transport_ * umol_to_mol * (atm_kpa_*kPa_to_Pa) * stom_cond_CO2_ * (-ca_ + c2 + 2 * gamma_ * umol_per_mol_to_Pa);
    double second_term = 16 * pow(stom_cond_CO2_, 2);
    double third_term = pow((ca_ + c2),2);
    double fourth_term = pow(electron_transport_, 2) * pow(umol_to_mol,2) * pow(atm_kpa_*kPa_to_Pa, 2);
    double fifth_term = 4*ca_*stom_cond_CO2_;
    double sixth_term = 4*c2*stom_cond_CO2_;
    double seventh_term = electron_transport_*umol_to_mol*(atm_kpa_*kPa_to_Pa);
    double eigth_term = 8*stom_cond_CO2_;

    return ci_ = (sqrt(first_term + second_term*third_term+ fourth_term) + fifth_term - sixth_term - seventh_term)/eigth_term;    
}

void Leaf::set_leaf_states_rates_from_psi_stem_analytical(double psi_stem) {

if (psi_soil_ >= psi_stem){
    stom_cond_CO2_ = 0;
    ci_ = gamma_*umol_per_mol_to_Pa;
      } else{
        psi_stem_to_ci_analytical(psi_stem);
    }

    assim_colimited_ = assim_colimited_analytical(ci_);    

    transpiration_ = stom_cond_CO2_ * H2O_CO2_stom_diff_ratio * atm_vpd / kg_to_mol_h2o / atm_kpa_;   
}

// Hydraulic cost equations

// Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::hydraulic_cost_Sperry(double psi_stem) {
  double k_l_soil_ = leaf_specific_conductance_max_ * proportion_of_conductivity(psi_soil_);
  double k_l_stem_ = leaf_specific_conductance_max_ * proportion_of_conductivity(psi_stem);
  
  hydraulic_cost_ = k_l_soil_ - k_l_stem_;
  
  return hydraulic_cost_;
}

double Leaf::hydraulic_cost_Bartlett(double psi_stem) {

hydraulic_cost_ = beta1 * sapwood_volume_per_leaf_area_ * pow((1 - proportion_of_conductivity(psi_stem)), beta2);

return hydraulic_cost_;
}

double Leaf::hydraulic_cost_TF(double psi_stem) {
hydraulic_cost_ = 1e6 * 
    hk_s /(365*24*60*60)* 
    (1/a_bio_) * 
    rho_ * sapwood_volume_per_leaf_area_ * pow((1 - proportion_of_conductivity(psi_stem)), beta2);

return hydraulic_cost_;
}

// Profit functions

double Leaf::profit_psi_stem_Sperry(double psi_stem) {

set_leaf_states_rates_from_psi_stem(psi_stem);

  double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_ * hydraulic_cost_;
}

double Leaf::profit_psi_stem_Bartlett(double psi_stem) {

set_leaf_states_rates_from_psi_stem(psi_stem);

  double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_Bartlett(psi_stem);

  return benefit_ - hydraulic_cost_;
}

double Leaf::profit_psi_stem_TF(double psi_stem) {
set_leaf_states_rates_from_psi_stem(psi_stem);

  double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_TF(psi_stem);

  return benefit_ - hydraulic_cost_;
}

double Leaf::profit_psi_stem_Bartlett_analytical(double psi_stem) {

set_leaf_states_rates_from_psi_stem_analytical(psi_stem);
  double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_Bartlett(psi_stem);

  return benefit_ - hydraulic_cost_;
}


double Leaf::profit_psi_stem_Sperry_analytical(double psi_stem) {

  set_leaf_states_rates_from_psi_stem_analytical(psi_stem);

  double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_analytical_ * hydraulic_cost_;
}

double Leaf::profit_Sperry_ci(double ci_) {                                  
  double benefit_ =
      assim_colimited(ci_);
  double stom_cond_CO2_ = (benefit_ * umol_to_mol * atm_kpa_ * kPa_to_Pa)/(ca_ - ci_); 
  double transpiration_ = stom_cond_CO2_ * H2O_CO2_stom_diff_ratio * atm_vpd_ / kg_to_mol_h2o / atm_kpa_;
  
  double psi_stem = transpiration_to_psi_stem(transpiration_);
  double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_*hydraulic_cost_;
}

double Leaf::profit_Sperry_ci_analytical(double ci_) {                                  
  double benefit_ =
      assim_colimited_analytical(ci_);

  double stom_cond_CO2_ = (benefit_ * umol_to_mol * atm_kpa_ * kPa_to_Pa)/(ca_ - ci_); 
  transpiration_ = stom_cond_CO2_ * H2O_CO2_stom_diff_ratio * atm_vpd_ / kg_to_mol_h2o / atm_kpa_;  

  double psi_stem = transpiration_to_psi_stem(transpiration_);
  double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_analytical_*hydraulic_cost_;
}


//optimisation functions


// need docs on Golden Section Search.
void Leaf::optimise_psi_stem_Sperry() {

  double gr = (sqrt(5) + 1) / 2;
  opt_psi_stem_ = psi_soil_;


  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    return;
  }

  // optimise for stem water potential
    double bound_a = psi_soil_;
    double bound_b = psi_crit;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > GSS_tol_abs) {

      double profit_at_c =
          profit_psi_stem_Sperry(bound_c);

      double profit_at_d =
          profit_psi_stem_Sperry(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem_ = ((bound_b + bound_a) / 2);
    profit_ = profit_psi_stem_Sperry(opt_psi_stem_);

  }
  
void Leaf::optimise_psi_stem_Sperry_analytical() {

  double gr = (sqrt(5) + 1) / 2;
  opt_psi_stem_ = psi_soil_;


  if (psi_soil_ > psi_crit){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;    
    return;  
    }

  // optimise for stem water potential
    double bound_a = psi_soil_;
    double bound_b = psi_crit;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > GSS_tol_abs) {

      double profit_at_c =
          profit_psi_stem_Sperry_analytical(bound_c);

      double profit_at_d =
          profit_psi_stem_Sperry_analytical(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }


    opt_psi_stem_ = ((bound_b + bound_a) / 2);
    profit_ = profit_psi_stem_Sperry_analytical(opt_psi_stem_);

  }

void Leaf::optimise_ci_Sperry(double max_ci) {

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_ci_ = gamma_*umol_per_mol_to_Pa;
    profit_ = 0.0;
    transpiration_ = 0.0;
    return;
  }

  double gr = (sqrt(5) + 1) / 2;

  
  // optimise for stem water potential
    double bound_a = gamma_*umol_per_mol_to_Pa;
    double bound_b = max_ci;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;
    while (abs(bound_b - bound_a) > GSS_tol_abs) {      

      double profit_at_c = profit_Sperry_ci(bound_c);

      double profit_at_d = profit_Sperry_ci(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_ci_ = ((bound_b + bound_a) / 2);
    profit_ = profit_Sperry_ci(opt_ci_);

  
    return;

  }

void Leaf::optimise_psi_stem_Bartlett() {

  double gr = (sqrt(5) + 1) / 2;
  opt_psi_stem_ = psi_soil_;

  if (psi_soil_ > psi_crit){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    return;
  }

  // optimise for stem water potential
    double bound_a = psi_soil_;
    double bound_b = psi_crit;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > GSS_tol_abs) {

      double profit_at_c =
          profit_psi_stem_Bartlett(bound_c);

      double profit_at_d =
          profit_psi_stem_Bartlett(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem_ = ((bound_b + bound_a) / 2);
    profit_ = profit_psi_stem_Bartlett(opt_psi_stem_);

    return;
      // }
  }

void Leaf::optimise_psi_stem_TF() {

  double gr = (sqrt(5) + 1) / 2;
  opt_psi_stem_ = psi_soil_;

  if (psi_soil_ > psi_crit){
    profit_ = profit_psi_stem_TF(psi_soil_);
    return;
  }

  // optimise for stem water potential
    double bound_a = psi_soil_;
    double bound_b = psi_crit;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;
    while (abs(bound_b - bound_a) > GSS_tol_abs) {

      double profit_at_c =
          profit_psi_stem_TF(bound_c);

      double profit_at_d =
          profit_psi_stem_TF(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem_ = ((bound_b + bound_a) / 2);
    profit_ = profit_psi_stem_TF(opt_psi_stem_);

    return;
  }

} // namespace plant
