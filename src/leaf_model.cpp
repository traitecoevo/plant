#include <plant/leaf_model.h>
#include <cmath>

namespace plant {
Leaf::Leaf(double vcmax, double c, double b,
           double psi_crit, // derived from b and c
           double epsilon_leaf, double beta1, double beta2, double jmax, double hydraulic_turnover)
    : vcmax(vcmax), // umol m^-2 s^-1 
    c(c), //unitless
    b(b), //-MPa
    psi_crit(psi_crit), //-MPa 
    beta1(beta1),
    beta2(beta2),
    epsilon_leaf(epsilon_leaf), //tolerance value 
    jmax(jmax),
    hydraulic_turnover(hydraulic_turnover), // s^-1
    leaf_specific_conductance_max_(NA_REAL), //kg m^-2 s^-1 MPa^-1 
    sapwood_volume_per_leaf_area_ (NA_REAL),
    PPFD_(NA_REAL), //? 
    atm_vpd_(NA_REAL), //kPa 
    ca_(NA_REAL), //Pa
    k_s_(NA_REAL), // yr ^ -1
    rho_(NA_REAL), //kg m^-3
    a_bio_(NA_REAL), //kg mol^-1
    psi_soil_(NA_REAL), //-MPa 
    ci_(NA_REAL), // Pa
    stom_cond_CO2_(NA_REAL), //mol Co2 m^-2 s^-1 
    electron_transport_(NA_REAL), //electron transport rate 
    assim_colimited_(NA_REAL), // umol C m^-2 s^-1 
    transpiration_(NA_REAL), // kg m^-2 s^-1 
    lambda_(NA_REAL), // umol C m^-2 s^-1 kg^-1 m^2 s^1
    lambda_analytical_(NA_REAL), // umol C m^-2 s^-1 kg^-1 m^2 s^1
    hydraulic_cost_(NA_REAL), // umol C m^-2 s^-1 
    profit_(NA_REAL), // umol C m^-2 s^-1 
    opt_psi_stem_(NA_REAL), //-MPa 
    opt_ci_(NA_REAL), //Pa 
    count(NA_REAL), 
    GSS_count(NA_REAL) {
      setup_transpiration(100);
}

//sets various parameters which are constant for a given node at a given time

void Leaf::set_physiology(double k_s, double rho, double a_bio, double PPFD, double psi_soil, double leaf_specific_conductance_max, double atm_vpd, double ca, double sapwood_volume_per_leaf_area) {
   k_s_ = k_s;
   rho_ = rho;
   a_bio_ = a_bio;
   atm_vpd_ = atm_vpd;
   PPFD_ = PPFD;
   psi_soil_ = psi_soil;
   leaf_specific_conductance_max_ = leaf_specific_conductance_max;
   sapwood_volume_per_leaf_area_ = sapwood_volume_per_leaf_area;
   ca_ = ca;
   electron_transport_ = electron_transport();

// set lambda, if psi_soil is higher than psi_crit, then set to 0. Currently doing both the numerical and analytical version. Ideally would do one.
  if(psi_soil >= psi_crit){
    lambda_ = 0;
    lambda_analytical_ = 0;

  } else {
    set_leaf_states_rates_from_psi_stem(psi_crit);
    lambda_ = assim_colimited(ci_) / hydraulic_cost_Sperry(psi_crit);
    set_leaf_states_rates_from_psi_stem_analytical(psi_crit);
    lambda_analytical_ = assim_colimited_analytical(ci_) / hydraulic_cost_Sperry(psi_crit);
  }
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
  double step = (psi_crit+psi_crit*0.1)/resolution;
  
  for (double psi_spline = 0.0 + step; psi_spline <= (psi_crit+psi_crit*0.1); psi_spline += step) {
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
  // std::cout << "E_stem" << transpiration_from_psi.eval(psi_stem) << std::endl;
  
  return leaf_specific_conductance_max_ * (transpiration_from_psi.eval(psi_stem) - transpiration_from_psi.eval(psi_soil_));

  
}

// converts a known transpiration to its corresponding psi_stem, returns -MPa
double Leaf::transpiration_to_psi_stem(double transpiration_) {
  // integration of proportion_of_conductivity over [psi_soil_, psi_stem]
  double E_psi_stem = transpiration_/leaf_specific_conductance_max_ +  transpiration_from_psi.eval(psi_soil_);

  return psi_from_transpiration.eval(E_psi_stem);
  }

// returns stomatal conductance to CO2, mol C m^-2 LA s^-1
double Leaf::stom_cond_CO2(double psi_stem) {
  double transpiration_ = transpiration(psi_stem);
  return atm_kpa * transpiration_ * kg_to_mol_h2o / atm_vpd_ / 1.6;
}


// biochemical photosynthesis model equations
//ensure that units of PPFD_ actually correspond to something real.
// electron trnansport rate based on light availability and vcmax assuming co-limitation hypothesis
double Leaf::electron_transport() {
  double electron_transport_ = (a * PPFD_ + jmax - sqrt(pow(a * PPFD_ + jmax, 2) - 4 * curv_fact * a * PPFD_ * jmax)) / (2 * curv_fact); // check brackets are correct
  return electron_transport_;           
}

//calculate the rubisco-limited assimilation rate, returns umol m^-2 s^-1
double Leaf::assim_rubisco_limited(double ci_) {
  return (vcmax * (ci_ - gamma_25 * umol_per_mol_to_Pa)) / (ci_ + km_25);
}

//calculate the light-limited assimilation rate, returns umol m^-2 s^-1
double Leaf::assim_electron_limited(double ci_) {
  
  return electron_transport_ / 4 *
  ((ci_ - gamma_25 * umol_per_mol_to_Pa) / (ci_ + 2 * gamma_25 * umol_per_mol_to_Pa));
}

// returns co-limited assimilation umol m^-2 s^-1
double Leaf::assim_colimited(double ci_) {
  
  double assim_rubisco_limited_ = assim_rubisco_limited(ci_);
  double assim_electron_limited_ = assim_electron_limited(ci_);
  
  // double R_d = vcmax * 0.015;
  // no dark respiration included at the moment
  return (assim_rubisco_limited_ + assim_electron_limited_ - sqrt(pow(assim_rubisco_limited_ + assim_electron_limited_, 2) - 4 * 0.98 * assim_rubisco_limited_ * assim_electron_limited_)) /
             (2 * 0.98);
}

// returns co-limited assimilation based only on light-limited transport rate and empirically-parameterised smothing term, returns umol m^-2 s^-1
double Leaf::assim_colimited_analytical(double ci_) {

  double c2 = 13.13652;

  // no dark respiration included at the moment
  
  return electron_transport_ / 4 *
         ((ci_ - gamma_25 * umol_per_mol_to_Pa) /
          (ci_ + c2));
}

// A - gc curves

// returns difference between co-limited assimilation and stom_cond_CO2, to be minimised (umol m^-2 s^-1)
double Leaf::assim_minus_stom_cond_CO2(double x, double psi_stem) {

  double assim_colimited_x_ = assim_colimited(x);

  double stom_cond_CO2_x_ = stom_cond_CO2(psi_stem);

  return assim_colimited_x_ * umol_to_mol -
         (stom_cond_CO2_x_ * (ca_ - x) / (atm_kpa * kPa_to_Pa));
}

// converts psi stem to ci, used to find ci which makes A(ci) = gc(ca - ci)
double Leaf::psi_stem_to_ci(double psi_stem) {
  // not clear what x is here
  
  auto target = [&](double x) mutable -> double {
    return assim_minus_stom_cond_CO2(x, psi_stem);
  };

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  return ci_ = util::uniroot(target, 0, ca_, 1e-6, 1000);
}

// given psi_stem, find assimilation, transpiration and stomal conductance to c02
void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem) {
  
  if (psi_soil_ >= psi_stem){
    ci_ = gamma_25*umol_per_mol_to_Pa;
    stom_cond_CO2_ = 0;
    } else{
      ci_ = psi_stem_to_ci(psi_stem);
      }
  
  assim_colimited_ = assim_colimited(ci_);

  transpiration_ = transpiration(psi_stem);
  stom_cond_CO2_ = atm_kpa * transpiration_ * kg_to_mol_h2o / atm_vpd_ / 1.6;

}

//this set of equations can make values greater than 40 (ca) when PPFD is extremely low
double Leaf::psi_stem_to_ci_analytical(double psi_stem) {

    stom_cond_CO2_ = stom_cond_CO2(psi_stem);
    double c2 = 13.13652;
        
    double first_term = 8 * electron_transport_ * umol_to_mol * (atm_kpa*kPa_to_Pa) * stom_cond_CO2_ * (-ca_ + c2 + 2 * gamma_25 * umol_per_mol_to_Pa);
    double second_term = 16 * pow(stom_cond_CO2_, 2);
    double third_term = pow((ca_ + c2),2);
    double fourth_term = pow(electron_transport_, 2) * pow(umol_to_mol,2) * pow(atm_kpa*kPa_to_Pa, 2);
    double fifth_term = 4*ca_*stom_cond_CO2_;
    double sixth_term = 4*c2*stom_cond_CO2_;
    double seventh_term = electron_transport_*umol_to_mol*(atm_kpa*kPa_to_Pa);
    double eigth_term = 8*stom_cond_CO2_;

    return ci_ = (sqrt(first_term + second_term*third_term+ fourth_term) + fifth_term - sixth_term - seventh_term)/eigth_term;    
}

void Leaf::set_leaf_states_rates_from_psi_stem_analytical(double psi_stem) {

if (psi_soil_ >= psi_stem){
    stom_cond_CO2_ = 0;
    ci_ = gamma_25*umol_per_mol_to_Pa;
      } else{
        psi_stem_to_ci_analytical(psi_stem);
    }

    assim_colimited_ = assim_colimited_analytical(ci_);    

    transpiration_ = stom_cond_CO2_ * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;   
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

hydraulic_cost_ = beta1 * sapwood_volume_per_leaf_area_ * pow((proportion_of_conductivity(psi_soil_) - proportion_of_conductivity(psi_stem)), beta2);

return hydraulic_cost_;
}

double Leaf::hydraulic_cost_TF(double psi_stem) {

hydraulic_cost_ = 1e6 * k_s_ /(365*24*60*60) * (1/a_bio_) * rho_ * sapwood_volume_per_leaf_area_ * pow((proportion_of_conductivity(psi_soil_) - proportion_of_conductivity(psi_stem)), beta2);

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
  double stom_cond_CO2_ = (benefit_ * umol_to_mol * atm_kpa * kPa_to_Pa)/(ca_ - ci_); 
  double transpiration_ = stom_cond_CO2_ * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;
  
  double psi_stem = transpiration_to_psi_stem(transpiration_);
  double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_*hydraulic_cost_;
}

double Leaf::profit_Sperry_ci_analytical(double ci_) {                                  
  double benefit_ =
      assim_colimited_analytical(ci_);

  double stom_cond_CO2_ = (benefit_ * umol_to_mol * atm_kpa * kPa_to_Pa)/(ca_ - ci_); 
  transpiration_ = stom_cond_CO2_ * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;  

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

    double delta_crit = 1e-07;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;
GSS_count = 0;

    while (abs(bound_b - bound_a) > delta_crit) {
GSS_count +=1 ;

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

    double delta_crit = 1e-2;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

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

void Leaf::optimise_psi_stem_Sperry_Newton(double psi_guess) {

  count = 0;
  int max_psi;
  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_psi_stem_ = psi_soil_;
    profit_ = 0.0;
    transpiration_ = 0.0;
    stom_cond_CO2_ = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.01; 
 
  double psi_stem_initial;

  int finished=1;

  opt_psi_stem_ = psi_guess;


  if (R_IsNA(opt_psi_stem_) | (opt_psi_stem_ > psi_crit) | (opt_psi_stem_ < psi_soil_)){

  max_psi = 1;

  opt_psi_stem_ = psi_crit - diff_value;

  }

  while(finished == 1){

    count += 1;

    if((count > 15) | ((R_IsNA(opt_psi_stem_) | (opt_psi_stem_ > psi_crit) | (opt_psi_stem_ < psi_soil_)) && max_psi == 3)){

      optimise_psi_stem_Sperry();
      return;
    }

  if ((R_IsNA(opt_psi_stem_) | (opt_psi_stem_ > psi_crit) | (opt_psi_stem_ < psi_soil_))  && max_psi == 1){

    count = 0;
    max_psi = 2;

    opt_psi_stem_ = psi_soil_ + diff_value;

  }

    if ((R_IsNA(opt_psi_stem_) | (opt_psi_stem_ > psi_crit) | (opt_psi_stem_ < psi_soil_)) && max_psi == 2){


    count = 0;

    max_psi = 3;

    opt_psi_stem_ = (psi_soil_ + psi_crit)/2;

  }

    psi_stem_initial = opt_psi_stem_;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;
    
    y_2 = profit_psi_stem_Sperry(x_2);

    y_1 = profit_psi_stem_Sperry(x_1);

    y_0 = profit_psi_stem_Sperry(psi_stem_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    opt_psi_stem_ = psi_stem_initial -  first_dev/sec_dev;

    if(abs(opt_psi_stem_ - psi_stem_initial) < epsilon_leaf){
      finished = 0;
    }

    profit_ = y_0;
  }
    return;
  }


void Leaf::optimise_psi_stem_Sperry_Newton_analytical(double psi_guess) {

  count = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_psi_stem_ = psi_soil_;
    profit_ = 0.0;
    transpiration_ = 0.0;
    stom_cond_CO2_ = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.001; 
 
  double psi_stem_initial;

  int finished=1;

  opt_psi_stem_ = psi_guess;

  if (R_IsNA(opt_psi_stem_) | (opt_psi_stem_ > psi_crit)){

  opt_psi_stem_ = psi_crit - diff_value;

  }

  while(finished == 1){

    count += 1;

  if (R_IsNA(opt_psi_stem_) | (opt_psi_stem_ > psi_crit)){

  opt_psi_stem_ = psi_crit - diff_value;

  }

    if(count > 10){
      optimise_psi_stem_Sperry_analytical();

      // std::cout << "gss" << std::endl;
      return;
    }


    psi_stem_initial = opt_psi_stem_;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;
    
    y_2 = profit_psi_stem_Sperry_analytical(x_2);

    y_1 = profit_psi_stem_Sperry_analytical(x_1);

    y_0 = profit_psi_stem_Sperry_analytical(psi_stem_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    opt_psi_stem_ = psi_stem_initial -  first_dev/sec_dev;

    if(y_0 == 0.0 & y_1 == 0.0 & y_2 == 0.0){
    opt_psi_stem_ = psi_soil_;
    profit_ = 0.0;
    transpiration_ = 0.0;
    std::cout << "warning2" << std::endl;
    return;
    }


    if(abs(opt_psi_stem_ - psi_stem_initial) < epsilon_leaf){
      finished = 0;
    }

    profit_ = y_0;
  }
    return;
  }



void Leaf::optimise_ci_Sperry_analytical(double max_ci) {

  if (psi_soil_ > psi_crit){

    opt_ci_ = gamma_25*umol_per_mol_to_Pa;
    profit_ = 0.0;
    transpiration_ = 0.0;
    stom_cond_CO2_ = 0.0;
    return;
  }
  
  double gr = (sqrt(5) + 1) / 2;

  
  // optimise for stem water potential
    double bound_a = gamma_25*umol_per_mol_to_Pa;
    double bound_b = max_ci;

    double delta_crit = 0.01;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c = profit_Sperry_ci_analytical(bound_c);

      double profit_at_d = profit_Sperry_ci_analytical(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_ci_ = ((bound_b + bound_a) / 2);
    profit_ = profit_Sperry_ci_analytical(opt_ci_);

  
    return;

  }



void Leaf::optimise_ci_Sperry(double max_ci) {

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_ci_ = gamma_25*umol_per_mol_to_Pa;
    profit_ = 0.0;
    transpiration_ = 0.0;
    return;
  }

  double gr = (sqrt(5) + 1) / 2;

  
  // optimise for stem water potential
    double bound_a = gamma_25*umol_per_mol_to_Pa;
    double bound_b = max_ci;

    double delta_crit = 0.01;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;
GSS_count = 0;
    while (abs(bound_b - bound_a) > delta_crit) {
GSS_count += 1;
      

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

void Leaf::optimise_ci_Sperry_Newton_analytical(double ci_guess) {

  int count = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_ci_ = gamma_25*umol_per_mol_to_Pa;
    profit_ = 0.0;
    transpiration_ = 0.0;
    stom_cond_CO2_ = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.001; 
  // TODO: rename as delta_ci
  
  
  double ci_initial;

  int unfinished=1; 

  opt_ci_ = ci_guess;

// add in the max ci
  if (R_IsNA(opt_ci_)){
  set_leaf_states_rates_from_psi_stem_analytical(psi_crit);

  opt_ci_ = ci_ - diff_value;

  }

  while(unfinished == 1){

    count += 1;

    if(count > 5){
      set_leaf_states_rates_from_psi_stem_analytical(psi_crit);
      optimise_ci_Sperry_analytical(ci_);


      return;

    }

    ci_initial = opt_ci_;
    
    if (!util::is_finite(opt_ci_)) {
      util::stop("Detected NAN ci value");
    }

    x_1 = ci_initial - diff_value;
    x_2 = ci_initial + diff_value;

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream
    double benefit_ = assim_colimited_analytical(x_2);


    double stom_cond_CO2_ = (benefit_ * umol_to_mol * atm_kpa * kPa_to_Pa)/(ca_ - x_2); 
    

    transpiration_ = stom_cond_CO2_ * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;  

//move to top

    double E_max = transpiration(psi_crit);

    if(((E_max < transpiration_) & (abs(transpiration_ - E_max) > 1e-15)) | (transpiration_ < 0)){

    set_leaf_states_rates_from_psi_stem_analytical(psi_crit); 

    opt_ci_ = ci_ - diff_value;


    continue;

    } else{

    double psi_stem = transpiration_to_psi_stem(transpiration_);

    double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem);
    y_2 = benefit_ - lambda_analytical_*hydraulic_cost_;

    }

    y_1 = profit_Sperry_ci_analytical(x_1);

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream

    y_0 = profit_Sperry_ci_analytical(ci_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

//added this in to account for situations where profit curve is flat and 0, likely due to no light
    if(y_0 == 0.0 & first_dev == 0.0 & sec_dev == 0.0){
    opt_ci_ = gamma_25*umol_per_mol_to_Pa;
    profit_ = 0.0;
    transpiration_ = 0.0;

    return;
    }

    opt_ci_ = ci_initial -  first_dev/sec_dev;

    if(abs(opt_ci_ - ci_initial) < (epsilon_leaf)){
     
      unfinished = 0;
    }

    profit_ = y_0;
  }

    return;
  }



void Leaf::optimise_ci_Sperry_Newton(double ci_guess) {

  count = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_ci_ = gamma_25*umol_per_mol_to_Pa;
    profit_ = 0.0;
    transpiration_ = 0.0;
    stom_cond_CO2_ = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.001; 
  // TODO: rename as delta_ci
  
  
  double ci_initial;

  int unfinished=1; 

  opt_ci_ = ci_guess;

// add in the max ci
  if (R_IsNA(opt_ci_)){
  set_leaf_states_rates_from_psi_stem(psi_crit);

  opt_ci_ = ci_ - diff_value;

  }

  while(unfinished == 1){

    count += 1;

    if(count > 10){

      set_leaf_states_rates_from_psi_stem(psi_crit);
      optimise_ci_Sperry(ci_);

      return;

    }

    ci_initial = opt_ci_;
    
    if (!util::is_finite(opt_ci_)) {
      util::stop("Detected NAN ci value");
    }

    x_1 = ci_initial - diff_value;
    x_2 = ci_initial + diff_value;

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream
    double benefit_ = assim_colimited(x_2);


    double stom_cond_CO2_ = (benefit_ * umol_to_mol * atm_kpa * kPa_to_Pa)/(ca_ - x_2); 
    

    transpiration_ = stom_cond_CO2_ * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;  

//move to top

    double E_max = transpiration(psi_crit);

    if(((E_max < transpiration_) & (abs(transpiration_ - E_max) > 1e-15)) | (transpiration_ < 0)){

    set_leaf_states_rates_from_psi_stem(psi_crit); 

    opt_ci_ = ci_ - diff_value;
    std::cout << "warning" << std::endl;


    continue;

    } 


    double psi_stem = transpiration_to_psi_stem(transpiration_);

    double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem);
    
    // start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream

    y_2 = benefit_ - lambda_*hydraulic_cost_;

    y_1 = profit_Sperry_ci(x_1);

    y_0 = profit_Sperry_ci(ci_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

//added this in to account for situations where profit curve is flat and 0, likely due to no light
    if(y_0 == 0.0 & first_dev == 0.0 & sec_dev == 0.0){
    opt_ci_ = gamma_25*umol_per_mol_to_Pa;
    profit_ = 0.0;
    transpiration_ = 0.0;
std::cout << "warning2" << std::endl;
    return;
    }

    opt_ci_ = ci_initial -  first_dev/sec_dev;


    if(abs(opt_ci_ - ci_initial) < (epsilon_leaf)){
     
      unfinished = 0;
    }

    profit_ = y_0;
  }

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

    double delta_crit = 1e-07;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;


GSS_count = 0;

    while (abs(bound_b - bound_a) > delta_crit) {
GSS_count +=1 ;

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

    if(profit_ < 0){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    opt_psi_stem_ = psi_soil_;
    return;
      }
  }



void Leaf::optimise_psi_stem_Bartlett_analytical() {
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

    double delta_crit = 1e-07;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;
GSS_count = 0;

    while (abs(bound_b - bound_a) > delta_crit) {
GSS_count +=1 ;

      double profit_at_c =
          profit_psi_stem_Bartlett_analytical(bound_c);

      double profit_at_d =
          profit_psi_stem_Bartlett_analytical(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem_ = ((bound_b + bound_a) / 2);
    profit_ = profit_psi_stem_Bartlett_analytical(opt_psi_stem_);
  }


void Leaf::optimise_psi_stem_TF() {

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

    double delta_crit = 1e-07;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;


GSS_count = 0;

    while (abs(bound_b - bound_a) > delta_crit) {
GSS_count +=1 ;

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

    if(profit_ < 0){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    opt_psi_stem_ = psi_soil_;
    return;
      }
  }

} // namespace plant