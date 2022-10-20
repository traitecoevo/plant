#include <plant/leaf_model.h>
#include <cmath>

namespace plant {
Leaf::Leaf(double vcmax, double p_50, double c, double b,
           double psi_crit, // derived from b and c
           double huber_value, double K_s, double epsilon_leaf)
    : vcmax(vcmax), p_50(p_50), c(c), b(b), psi_crit(psi_crit), huber_value(huber_value), K_s(K_s), epsilon_leaf(epsilon_leaf), 
      ci(NA_REAL), g_c(NA_REAL), A_lim_(NA_REAL), E(NA_REAL), profit(NA_REAL), j_(NA_REAL), lambda_(NA_REAL),
      lambda_analytical_(NA_REAL), PPFD_(NA_REAL), atm_vpd_(NA_REAL), ca_(NA_REAL),psi_soil_(NA_REAL), k_l_max_(NA_REAL), opt_psi_stem(NA_REAL), opt_ci(NA_REAL), method(NA_REAL),
      count(NA_REAL), GSS_count(NA_REAL) {
        setup_E_supply(100);
}

void Leaf::set_physiology(double PPFD, double psi_soil, double k_l_max, double atm_vpd, double ca) {
   atm_vpd_ = atm_vpd;
   PPFD_ = PPFD;
   psi_soil_ = psi_soil;
   k_l_max_ = k_l_max;
   ca_ = ca;
   j_ = calc_j();

  if(psi_soil >= psi_crit){
    lambda_ = 0;
    lambda_analytical_ = 0;

  } else {
    get_leaf_states_rates_from_psi_stem(psi_crit);
    lambda_ = A_lim(ci) / calc_hydraulic_cost_Sperry(psi_crit);
    get_leaf_states_rates_from_psi_stem_analytical(psi_crit);
    lambda_analytical_ = A_lim_analytical(ci) / calc_hydraulic_cost_Sperry(psi_crit);
  }
}

// transpiration supply functions

// integrates, returns proportion of conductance taken from hydraulic vulnerability curve (unitless)
double Leaf::calc_cond_vuln(double psi) const {
  return exp(-pow((psi / b), c));
}
// REMOVED k_l_max_
void Leaf::setup_E_supply(double resolution) {
  // integrate and accumulate results
  auto x_psi = std::vector<double>{0.0};  // {0.0}
  auto y_cumulative_E = std::vector<double>{0.0}; // {0.0}
  double step = (psi_crit+psi_crit*0.1)/resolution;
  
  for (double psi_spline = 0.0 + step; psi_spline <= (psi_crit+psi_crit*0.1); psi_spline += step) {
    double E_psi = step * ((calc_cond_vuln(psi_spline-step) + calc_cond_vuln(psi_spline))/2) + y_cumulative_E.back();
    x_psi.push_back(psi_spline); // x values for spline
    y_cumulative_E.push_back(E_psi); // y values for spline

}
// setup interpolator
E_from_psi.init(x_psi, y_cumulative_E);
E_from_psi.set_extrapolate(false);

psi_from_E.init(y_cumulative_E, x_psi);
psi_from_E.set_extrapolate(false);
}

// replace f with some other function, returns E kg m^-2 s^-1

double Leaf::calc_E_supply_full_integration(double psi_stem) {
  std::function<double(double)> f;
  f = [&](double psi) -> double { return calc_cond_vuln(psi); };
  
  return k_l_max_ * integrator.integrate(f, psi_soil_, psi_stem);
 }

//returns kg h20 s^-1 m^-2 LA
double Leaf::calc_E_supply(double psi_stem) {
  // integration of calc_cond_vuln over [psi_soil_, psi_stem]
  // std::cout << "E_stem" << E_from_psi.eval(psi_stem) << std::endl;
  
  return k_l_max_ * (E_from_psi.eval(psi_stem) - E_from_psi.eval(psi_soil_));

  
}

double Leaf::convert_E_from_ci_to_psi_stem(double E_ci) {
  // integration of calc_cond_vuln over [psi_soil_, psi_stem]
  double E_psi_stem = E_ci/k_l_max_ +  E_from_psi.eval(psi_soil_);

  return psi_from_E.eval(E_psi_stem);
  }

// returns g_c, mol C m^-2 LA s^-1
double Leaf::calc_g_c(double psi_stem) {
  double E_supply = calc_E_supply(psi_stem);
  return atm_kpa * E_supply * kg_to_mol_h2o / atm_vpd_ / 1.6;
}


// biochemical photosynthesis model equations
//ensure that units of PPFD_ actually correspond to something real.
double Leaf::calc_j() {
  double jmax = vcmax * vcmax_25_to_jmax_25;
  double j = (a * PPFD_ + jmax - sqrt(pow(a * PPFD_ + jmax, 2) - 4 * curv_fact * a * PPFD_ * jmax)) / (2 * curv_fact); // check brackets are correct

  return j;           
}

// returns A_c umol m^-2 s^-1
double Leaf::calc_A_c(double ci_) {
  return (vcmax * (ci_ - gamma_25 * umol_per_mol_to_Pa)) / (ci_ + km_25);
}

// returns A_j umol m^-2 s^-1
double Leaf::calc_A_j(double ci_) {
  
  return j_ / 4 *
  ((ci_ - gamma_25 * umol_per_mol_to_Pa) / (ci_ + 2 * gamma_25 * umol_per_mol_to_Pa));
}

// returns co-limited assimilation umol m^-2 s^-1
double Leaf::A_lim(double ci_) {
  
  double A_c = calc_A_c(ci_);
  double A_j = calc_A_j(ci_);
  
  // double R_d = vcmax * 0.015;
  // no dark respiration included at the moment
  return (A_c + A_j - sqrt(pow(A_c + A_j, 2) - 4 * 0.98 * A_c * A_j)) /
             (2 * 0.98);
}

// returns co-limited assimilation umol m^-2 s^-1
double Leaf::A_lim_analytical(double ci_) {

  double c2 = 13.13652;

  // no dark respiration included at the moment
  // double j = calc_j();
  
  return j_ / 4 *
         ((ci_ - gamma_25 * umol_per_mol_to_Pa) /
          (ci_ + c2));
}

// A - gc curves

// returns difference between co-limited assimilation and g_c, to be minimised (umol m^-2 s^-1)
double Leaf::diff_ci(double x, double psi_stem) {

  double A_lim_x = A_lim(x);

  double g_c_ = calc_g_c(psi_stem);

  return A_lim_x * umol_to_mol -
         (g_c_ * (ca_ - x) / (atm_kpa * kPa_to_Pa));
}


double Leaf::convert_psi_stem_to_ci(double psi_stem) {
  // not clear what x is here
  
  auto target = [&](double x) mutable -> double {
    return diff_ci(x, psi_stem);
  };

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  return ci = util::uniroot(target, 0, ca_, 1e-6, 1000);
}

void Leaf::get_leaf_states_rates_from_psi_stem(double psi_stem) {
  
  if (psi_soil_ >= psi_stem){
    ci = gamma_25*umol_per_mol_to_Pa;
    g_c = 0;
    } else{
      ci = convert_psi_stem_to_ci(psi_stem);
      }
  
  A_lim_ = A_lim(ci);
  // g_c = calc_g_c(psi_stem);
  // E = g_c * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;

  
  E = calc_E_supply(psi_stem);
  g_c = atm_kpa * E * kg_to_mol_h2o / atm_vpd_ / 1.6;

}

//this set of equations can make values greater than 40 (ca) when PPFD is extremely low
double Leaf::convert_psi_stem_to_ci_analytical(double psi_stem) {

    g_c = calc_g_c(psi_stem);
    double c2 = 13.13652;
        
    double first_term = 8 * j_ * umol_to_mol * (atm_kpa*kPa_to_Pa) * g_c * (-ca_ + c2 + 2 * gamma_25 * umol_per_mol_to_Pa);
    double second_term = 16 * pow(g_c, 2);
    double third_term = pow((ca_ + c2),2);
    double fourth_term = pow(j_, 2) * pow(umol_to_mol,2) * pow(atm_kpa*kPa_to_Pa, 2);
    double fifth_term = 4*ca_*g_c;
    double sixth_term = 4*c2*g_c;
    double seventh_term = j_*umol_to_mol*(atm_kpa*kPa_to_Pa);
    double eigth_term = 8*g_c;

    return ci = (sqrt(first_term + second_term*third_term+ fourth_term) + fifth_term - sixth_term - seventh_term)/eigth_term;    
}

void Leaf::get_leaf_states_rates_from_psi_stem_analytical(double psi_stem) {

if (psi_soil_ >= psi_stem){
    g_c = 0;
    ci = gamma_25*umol_per_mol_to_Pa;
      } else{
        convert_psi_stem_to_ci_analytical(psi_stem);
    }

    A_lim_ = A_lim_analytical(ci);    
    // g_c = calc_g_c(psi_stem);

    E = g_c * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;   
}

// Hydraulic cost equations

// Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::calc_hydraulic_cost_Sperry(double psi_stem) {
  double k_l_soil_ = k_l_max_ * calc_cond_vuln(psi_soil_);
  double k_l_stem_ = k_l_max_ * calc_cond_vuln(psi_stem);

  return k_l_soil_ - k_l_stem_;
}

// Profit functions

double Leaf::profit_psi_stem_Sperry(double psi_stem) {

get_leaf_states_rates_from_psi_stem(psi_stem);

  double benefit_ = A_lim_;
  double cost_ = calc_hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_ * cost_;
}

double Leaf::profit_psi_stem_Sperry_analytical(double psi_stem) {

  get_leaf_states_rates_from_psi_stem_analytical(psi_stem);

  double benefit_ = A_lim_;
  double cost_ = calc_hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_analytical_ * cost_;
}

double Leaf::calc_profit_Sperry_ci(double c_i) {                                  
  double benefit_ =
      A_lim(c_i);
  double g_c_ci = (benefit_ * umol_to_mol * atm_kpa * kPa_to_Pa)/(ca_ - c_i); 
  double E_ci = g_c_ci * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;
  
  double psi_stem = convert_E_from_ci_to_psi_stem(E_ci);
  double cost_ = calc_hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_*cost_;
}

double Leaf::calc_profit_Sperry_ci_analytical(double c_i) {                                  
  double benefit_ =
      A_lim_analytical(c_i);

  double g_c_ci = (benefit_ * umol_to_mol * atm_kpa * kPa_to_Pa)/(ca_ - c_i); 
  E = g_c_ci * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;  

  double psi_stem = convert_E_from_ci_to_psi_stem(E);
  double cost_ = calc_hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_analytical_*cost_;
}


//optimisation functions


// need docs on Golden Section Search.
void Leaf::optimise_psi_stem_Sperry() {

  double gr = (sqrt(5) + 1) / 2;
  opt_psi_stem = psi_soil_;


  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    profit = 0;
    E = 0;
    g_c = 0;
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

    opt_psi_stem = ((bound_b + bound_a) / 2);
    profit = profit_psi_stem_Sperry(opt_psi_stem);

    // std::cout << "E" << E << std::endl;

  }
  
void Leaf::optimise_psi_stem_Sperry_analytical() {

  double gr = (sqrt(5) + 1) / 2;
  opt_psi_stem = psi_soil_;


  if (psi_soil_ > psi_crit){
    profit = 0;
    E = 0;
    g_c = 0;    
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
    opt_psi_stem = ((bound_b + bound_a) / 2);
    profit = profit_psi_stem_Sperry_analytical(opt_psi_stem);

  }



void Leaf::optimise_psi_stem_Sperry_Newton(double psi_guess) {

  count = 0;
  int max_psi;
  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_psi_stem = psi_soil_;
    profit = 0.0;
    E = 0.0;
    g_c = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.01; 
 
  double psi_stem_initial;

  int finished=1;

  opt_psi_stem = psi_guess;

  method = 0;

  if (R_IsNA(opt_psi_stem) | (opt_psi_stem > psi_crit) | (opt_psi_stem < psi_soil_)){

  max_psi = 1;

  opt_psi_stem = psi_crit - diff_value;

  }

  while(finished == 1){

    count += 1;

    if((count > 15) | ((R_IsNA(opt_psi_stem) | (opt_psi_stem > psi_crit) | (opt_psi_stem < psi_soil_)) && max_psi == 3)){
      method = 1;

      optimise_psi_stem_Sperry();
      return;
    }

  if ((R_IsNA(opt_psi_stem) | (opt_psi_stem > psi_crit) | (opt_psi_stem < psi_soil_))  && max_psi == 1){

    count = 0;
    max_psi = 2;

    opt_psi_stem = psi_soil_ + diff_value;

  }

    if ((R_IsNA(opt_psi_stem) | (opt_psi_stem > psi_crit) | (opt_psi_stem < psi_soil_)) && max_psi == 2){


    count = 0;

    max_psi = 3;

    opt_psi_stem = (psi_soil_ + psi_crit)/2;

  }

    psi_stem_initial = opt_psi_stem;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;
    
    y_2 = profit_psi_stem_Sperry(x_2);

    y_1 = profit_psi_stem_Sperry(x_1);

    y_0 = profit_psi_stem_Sperry(psi_stem_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    opt_psi_stem = psi_stem_initial -  first_dev/sec_dev;

    if(abs(opt_psi_stem - psi_stem_initial) < epsilon_leaf){
      finished = 0;
    }

    profit = y_0;
  }
    return;
  }


void Leaf::optimise_psi_stem_Sperry_Newton_analytical(double psi_guess) {

  count = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_psi_stem = psi_soil_;
    profit = 0.0;
    E = 0.0;
    g_c = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.001; 
 
  double psi_stem_initial;

  int finished=1;

  opt_psi_stem = psi_guess;

  if (R_IsNA(opt_psi_stem) | (opt_psi_stem > psi_crit)){

  opt_psi_stem = psi_crit - diff_value;

  }

  while(finished == 1){

    count += 1;

  if (R_IsNA(opt_psi_stem) | (opt_psi_stem > psi_crit)){

  opt_psi_stem = psi_crit - diff_value;

  }

    if(count > 10){
      optimise_psi_stem_Sperry_analytical();

      // std::cout << "gss" << std::endl;
      return;
    }


    psi_stem_initial = opt_psi_stem;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;
    
    y_2 = profit_psi_stem_Sperry_analytical(x_2);

    y_1 = profit_psi_stem_Sperry_analytical(x_1);

    y_0 = profit_psi_stem_Sperry_analytical(psi_stem_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    opt_psi_stem = psi_stem_initial -  first_dev/sec_dev;

    if(y_0 == 0.0 & y_1 == 0.0 & y_2 == 0.0){
    opt_psi_stem = psi_soil_;
    profit = 0.0;
    E = 0.0;
    std::cout << "warning2" << std::endl;
    return;
    }


    if(abs(opt_psi_stem - psi_stem_initial) < epsilon_leaf){
      finished = 0;
    }

    profit = y_0;
  }
    return;
  }



void Leaf::optimise_ci_Sperry_analytical(double max_ci) {

  if (psi_soil_ > psi_crit){

    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0.0;
    E = 0.0;
    g_c = 0.0;
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

      double profit_at_c = calc_profit_Sperry_ci_analytical(bound_c);

      double profit_at_d = calc_profit_Sperry_ci_analytical(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_ci = ((bound_b + bound_a) / 2);
    profit = calc_profit_Sperry_ci_analytical(opt_ci);

  
    return;

  }



void Leaf::optimise_ci_Sperry(double max_ci) {

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0.0;
    E = 0.0;
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
      

      double profit_at_c = calc_profit_Sperry_ci(bound_c);

      double profit_at_d = calc_profit_Sperry_ci(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_ci = ((bound_b + bound_a) / 2);
    profit = calc_profit_Sperry_ci(opt_ci);

  
    return;

  }

void Leaf::optimise_ci_Sperry_Newton_analytical(double ci_guess) {

  int count = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0.0;
    E = 0.0;
    g_c = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.001; 
  // TODO: rename as delta_ci
  
  
  double ci_initial;

  int unfinished=1; 

  opt_ci = ci_guess;

// add in the max ci
  if (R_IsNA(opt_ci)){
  get_leaf_states_rates_from_psi_stem_analytical(psi_crit);

  opt_ci = ci - diff_value;

  }

  while(unfinished == 1){

    count += 1;

    if(count > 5){
      get_leaf_states_rates_from_psi_stem_analytical(psi_crit);
      optimise_ci_Sperry_analytical(ci);


      return;

    }

    ci_initial = opt_ci;
    
    if (!util::is_finite(opt_ci)) {
      util::stop("Detected NAN ci value");
    }

    x_1 = ci_initial - diff_value;
    x_2 = ci_initial + diff_value;

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream
    double benefit_ = A_lim_analytical(x_2);


    double g_c_ci = (benefit_ * umol_to_mol * atm_kpa * kPa_to_Pa)/(ca_ - x_2); 
    

    E = g_c_ci * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;  

//move to top

    double E_max = calc_E_supply(psi_crit);

    if(((E_max < E) & (abs(E - E_max) > 1e-15)) | (E < 0)){

    get_leaf_states_rates_from_psi_stem_analytical(psi_crit); 

    opt_ci = ci - diff_value;


    continue;

    } else{

    double psi_stem = convert_E_from_ci_to_psi_stem(E);

    double cost_ = calc_hydraulic_cost_Sperry(psi_stem);
    y_2 = benefit_ - lambda_analytical_*cost_;

    }

    y_1 = calc_profit_Sperry_ci_analytical(x_1);

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream

    y_0 = calc_profit_Sperry_ci_analytical(ci_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

//added this in to account for situations where profit curve is flat and 0, likely due to no light
    if(y_0 == 0.0 & first_dev == 0.0 & sec_dev == 0.0){
    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0.0;
    E = 0.0;

    return;
    }

    opt_ci = ci_initial -  first_dev/sec_dev;

    if(abs(opt_ci - ci_initial) < (epsilon_leaf)){
     
      unfinished = 0;
    }

    profit = y_0;
  }

    return;
  }



void Leaf::optimise_ci_Sperry_Newton(double ci_guess) {

  count = 0;

  method = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0.0;
    E = 0.0;
    g_c = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.001; 
  // TODO: rename as delta_ci
  
  
  double ci_initial;

  int unfinished=1; 

  opt_ci = ci_guess;

// add in the max ci
  if (R_IsNA(opt_ci)){
  get_leaf_states_rates_from_psi_stem(psi_crit);

  opt_ci = ci - diff_value;

  }

  while(unfinished == 1){

    count += 1;

    if(count > 10){

      method = 1;

      get_leaf_states_rates_from_psi_stem(psi_crit);
      optimise_ci_Sperry(ci);

      return;

    }

    ci_initial = opt_ci;
    
    if (!util::is_finite(opt_ci)) {
      util::stop("Detected NAN ci value");
    }

    x_1 = ci_initial - diff_value;
    x_2 = ci_initial + diff_value;

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream
    double benefit_ = A_lim(x_2);


    double g_c_ci = (benefit_ * umol_to_mol * atm_kpa * kPa_to_Pa)/(ca_ - x_2); 
    

    E = g_c_ci * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;  

//move to top

    double E_max = calc_E_supply(psi_crit);

    if(((E_max < E) & (abs(E - E_max) > 1e-15)) | (E < 0)){

    get_leaf_states_rates_from_psi_stem(psi_crit); 

    opt_ci = ci - diff_value;
    std::cout << "warning" << std::endl;


    continue;

    } 


    double psi_stem = convert_E_from_ci_to_psi_stem(E);

    double cost_ = calc_hydraulic_cost_Sperry(psi_stem);
    
    // start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream

    y_2 = benefit_ - lambda_*cost_;

    y_1 = calc_profit_Sperry_ci(x_1);

    y_0 = calc_profit_Sperry_ci(ci_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

//added this in to account for situations where profit curve is flat and 0, likely due to no light
    if(y_0 == 0.0 & first_dev == 0.0 & sec_dev == 0.0){
    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0.0;
    E = 0.0;
std::cout << "warning2" << std::endl;
    return;
    }

    opt_ci = ci_initial -  first_dev/sec_dev;


    if(abs(opt_ci - ci_initial) < (epsilon_leaf)){
     
      unfinished = 0;
    }

    profit = y_0;
  }

    return;
  }



} // namespace plant