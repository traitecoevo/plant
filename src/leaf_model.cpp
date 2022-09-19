#include <plant/leaf_model.h>
#include <cmath>

namespace plant {
Leaf::Leaf(double vcmax, double p_50, double c, double b,
           double psi_crit, // derived from b and c
           double beta, double beta_2, double huber_value, double K_s, double epsilon_leaf)
    : vcmax(vcmax), p_50(p_50), c(c), b(b), psi_crit(psi_crit), beta(beta),
      beta_2(beta_2), huber_value(huber_value), K_s(K_s), epsilon_leaf(epsilon_leaf), ci(NA_REAL), g_c(NA_REAL),
      A_lim(NA_REAL), E(NA_REAL), psi(NA_REAL), profit(NA_REAL), psi_stem_next(NA_REAL),j_(NA_REAL), c_i_next(NA_REAL), lambda_(NA_REAL),
      PPFD_(NA_REAL), atm_vpd_(NA_REAL), psi_soil_(NA_REAL), k_l_max_(NA_REAL), opt_psi_stem(NA_REAL), opt_ci(NA_REAL) {
    setup_E_supply(100);
    // setup_psi(100);

    // std::cout << "psi_crit" << psi_crit << std::endl;
}

void Leaf::set_physiology(double PPFD, double psi_soil, double k_l_max, double atm_vpd) {
   atm_vpd_ = atm_vpd;
   PPFD_ = PPFD;
   psi_soil_ = psi_soil;
   k_l_max_ = k_l_max;
   j_ = calc_j();

  if(psi_soil >= psi_crit){
    lambda_ = 0;
  } else {
    lambda_ = calc_assim_gross_one_line(psi_crit) / (k_l_max_ * calc_cond_vuln(psi_soil_) - k_l_max_ * calc_cond_vuln(psi_crit));
  }
  //  std::cout << "max_assim" << calc_assim_gross_one_line(psi_crit) << "first_term_denom" << k_l_max_ * calc_cond_vuln(psi_soil_) << "second_term_denom" << k_l_max_ * calc_cond_vuln(psi_crit) << std::endl;
}

// REMOVED k_l_max_
// integrates, returns conductivity at given psi_stem kg m^-2 s^-1 MPa^-1

double Leaf::calc_cond_vuln(double psi) const {
  // std::cout << "psi" << psi << "b" << b << "c" << c;
  return exp(-pow((psi / b), c));
}

// REMOVED k_l_max_
void Leaf::setup_E_supply(double resolution) {
    // integrate and accumulate results
    auto x_psi = std::vector<double>{0.0};  // {0.0}
    auto y_cumulative_E = std::vector<double>{0.0}; // {0.0}
    double step = (psi_crit+psi_crit*0.1)/resolution;
    for (double psi_spline = 0.0 + step; psi_spline <= (psi_crit+psi_crit*0.1); psi_spline += step) {
        // std::cout << "psi_spline" << psi_spline << std::endl;
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

double Leaf::calc_E_supply(double psi_stem) {
    // integration of calc_cond_vuln over [psi_soil_, psi_stem]
    return k_l_max_ * (E_from_psi.eval(psi_stem) - E_from_psi.eval(psi_soil_));
}

double Leaf::calc_psi_stem_ci(double E_ci) {
    // integration of calc_cond_vuln over [psi_soil_, psi_stem]
    double E_psi_stem = E_ci/k_l_max_ +  E_from_psi.eval(psi_soil_);
    return psi_from_E.eval(E_psi_stem);
    }

// double Leaf::calc_psi_from_E(double E_psi_stem) {
//     // integration of calc_cond_vuln over [psi_soil_, psi_stem]
//     return psi_from_E.eval(E_psi_stem);
//     }


// returns E kg m^-2 s^-1
double Leaf::calc_g_c(double psi_stem) {
  double E_supply = calc_E_supply(psi_stem);

  return atm_kpa * E_supply * kg_to_mol_h2o / atm_vpd_ / 1.6;
}

double Leaf::calc_j() {
  double jmax = vcmax * vcmax_25_to_jmax_25;
  double j = (a * PPFD_ + jmax -
              sqrt(pow(a * PPFD_ + jmax, 2) - 4 * curv_fact * a * PPFD_ * jmax)) /
             (2 * curv_fact); // check brackets are correct

  return j;           

}

// returns A_c umol m^-2 s^-1
double Leaf::calc_A_c(double ci_) {
  return (vcmax * (ci_ - gamma_25 * umol_per_mol_to_Pa)) / (ci_ + km_25);
}

// returns A_j umol m^-2 s^-1
double Leaf::calc_A_j(double ci_) {

  // double j = calc_j();
  return j_ / 4 *
         ((ci_ - gamma_25 * umol_per_mol_to_Pa) /
          (ci_ + 2 * gamma_25 * umol_per_mol_to_Pa));
}

// returns co-limited assimilation umol m^-2 s^-1
double Leaf::calc_A_lim(double ci_) {

  double A_c = calc_A_c(ci_);
  double A_j = calc_A_j(ci_);
  // double R_d = vcmax * 0.015;

  return (A_c + A_j - sqrt(pow(A_c + A_j, 2) - 4 * 0.98 * A_c * A_j)) /
             (2 * 0.98);
}

// returns co-limited assimilation umol m^-2 s^-1
double Leaf::calc_A_lim_one_line(double ci_) {

  double c2 = 13.13652;

  // double j = calc_j();
  return j_ / 4 *
         ((ci_ - gamma_25 * umol_per_mol_to_Pa) /
          (ci_ + c2));
}

// returns difference between co-limited assimilation and g_c, to be minimised (umol m^-2 s^-1)
double Leaf::diff_ci(double x, double psi_stem) {

  double A_lim_ = calc_A_lim(x);

  double g_c_ = calc_g_c(psi_stem);

  return A_lim_ * umol_per_mol_to_mol_per_mol -
         (g_c_ * (ca - x) / (atm_kpa * kPa_to_Pa));
}

// need to fill in tol and max_iteratiosn
double Leaf::calc_assim_gross(double psi_stem) {

 if (psi_soil_ == psi_stem){
    ci = gamma_25*umol_per_mol_to_Pa;
      } else{

  // not clear what x is here
  auto target = [&](double x) mutable -> double {
    return diff_ci(x, psi_stem);
  };

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  ci = util::uniroot(target, 0, ca, 1e-6, 1000);
      }
  A_lim = calc_A_lim(ci);

  g_c = calc_g_c(psi_stem);

  //E is in m^3 m-2 

  E = g_c * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;

  psi = psi_stem;

  return A_lim;
  
}

double Leaf::calc_assim_gross_one_line(double psi_stem) {

    if (psi_soil_ >= psi_stem){
    g_c = 0;
    ci = gamma_25*umol_per_mol_to_Pa;
      } else{
    g_c = calc_g_c(psi_stem);
    
    double c2 = 13.13652;
    
    // double j = calc_j();
    
    double first_term = 8 * j_ * umol_per_mol_to_mol_per_mol * (atm_kpa*kPa_to_Pa) * g_c * (-ca + c2 + 2 * gamma_25 * umol_per_mol_to_Pa);
    double second_term = 16 * pow(g_c, 2);
    double third_term = pow((ca + c2),2);
    double fourth_term = pow(j_, 2) * pow(umol_per_mol_to_mol_per_mol,2) * pow(atm_kpa*kPa_to_Pa, 2);
    double fifth_term = 4*ca*g_c;
    double sixth_term = 4*c2*g_c;
    double seventh_term = j_*umol_per_mol_to_mol_per_mol*(atm_kpa*kPa_to_Pa);
    double eigth_term = 8*g_c;

    ci = (sqrt(first_term + second_term*third_term+ fourth_term) + fifth_term - sixth_term - seventh_term)/eigth_term;
      } 
    A_lim = calc_A_lim_one_line(ci);
    
    E = g_c * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;
    
    psi = psi_stem;
    
    return A_lim;
}


// Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::calc_hydraulic_cost_Sperry(double psi_stem) {
  double k_l_soil_ = k_l_max_ * calc_cond_vuln(psi_soil_);
  double k_l_stem_ = k_l_max_ * calc_cond_vuln(psi_stem);

  return k_l_soil_ - k_l_stem_;
}

double Leaf::find_max_ci() {

  // not clear what x is he
  auto target = [&](double x) mutable -> double {
    return diff_ci(x, psi_crit);
  };

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  ci = util::uniroot(target, 0, ca, 1e-6, 1000);

  return ci;
  
}

//this set of equations can make values greater than 40 (ca) when PPFD is extremely low
double Leaf::find_max_ci_one_line() {
    
    g_c = calc_g_c(psi_crit);
    // std::cout << "psi_crit_find_max" << psi_crit << std::endl;
    double c2 = 13.13652;
    // double j = calc_j();
    
    double first_term = 8 * j_ * umol_per_mol_to_mol_per_mol * (atm_kpa*kPa_to_Pa) * g_c * (-ca + c2 + 2 * gamma_25 * umol_per_mol_to_Pa);
    double second_term = 16 * pow(g_c, 2);
    double third_term = pow((ca + c2),2);
    double fourth_term = pow(j_, 2) * pow(umol_per_mol_to_mol_per_mol,2) * pow(atm_kpa*kPa_to_Pa, 2);
    double fifth_term = 4*ca*g_c;
    double sixth_term = 4*c2*g_c;
    double seventh_term = j_*umol_per_mol_to_mol_per_mol*(atm_kpa*kPa_to_Pa);
    double eigth_term = 8*g_c;

    ci = (sqrt(first_term + second_term*third_term+ fourth_term) + fifth_term - sixth_term - seventh_term)/eigth_term;    
    return ci;
}

double Leaf::calc_profit_Sperry(double psi_stem) {

  double benefit_ =
      calc_assim_gross(psi_stem);
  double cost_ = calc_hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_ * cost_;
}

double Leaf::calc_profit_Sperry_one_line(double psi_stem) {

  double benefit_ =
      calc_assim_gross_one_line(psi_stem);
  double cost_ = calc_hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_ * cost_;
}

double Leaf::calc_profit_Sperry_ci(double c_i) {                                  
  double benefit_ =
      calc_A_lim(c_i);
  double g_c_ci = (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - c_i); 
  double E_ci = g_c_ci * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;
  // std::cout << "E_ci" << E_ci << "g_c" << g_c_ci << "c_i" << c_i << std::endl;
  
  double psi_stem = calc_psi_stem_ci(E_ci);
  double cost_ = calc_hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_*cost_;
}

double Leaf::calc_profit_Sperry_ci_one_line(double c_i) {                                  
  double benefit_ =
      calc_A_lim_one_line(c_i);

  double g_c_ci = (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - c_i); 

  E = g_c_ci * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;  

  double psi_stem = calc_psi_stem_ci(E);
  double cost_ = calc_hydraulic_cost_Sperry(psi_stem);
  return benefit_ - lambda_*cost_;
}


// need docs on Golden Section Search.
double Leaf::optimise_psi_stem_Sperry() {

  double gr = (sqrt(5) + 1) / 2;
  double opt_psi_stem = psi_soil_;


  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    return(opt_psi_stem);
  }

  // optimise for stem water potential
    double bound_a = psi_soil_;
    double bound_b = psi_crit;

    double delta_crit = 1e-3;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Sperry(bound_c);

      double profit_at_d =
          calc_profit_Sperry(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem = ((bound_b + bound_a) / 2);
  
    return opt_psi_stem;

  }
  
double Leaf::optimise_psi_stem_Sperry_one_line() {
    //  std::cout << "hello" << std::endl;

  double gr = (sqrt(5) + 1) / 2;
  double opt_psi_stem = psi_soil_;


  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    return(opt_psi_stem);
  }

  // optimise for stem water potential
    double bound_a = psi_soil_;
    double bound_b = psi_crit;

    double delta_crit = 0.01;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Sperry_one_line(bound_c);

      double profit_at_d =
          calc_profit_Sperry_one_line(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem = ((bound_b + bound_a) / 2);

  
    return opt_psi_stem;

  }

void Leaf::optimise_ci_Sperry_one_line(double max_ci) {
    //  std::cout << "hello" << std::endl;

  double gr = (sqrt(5) + 1) / 2;

  
  // optimise for stem water potential
    double bound_a = gamma_25*umol_per_mol_to_Pa;
    double bound_b = max_ci;

    double delta_crit = 0.01;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c = calc_profit_Sperry_ci_one_line(bound_c);

      double profit_at_d = calc_profit_Sperry_ci_one_line(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_ci = ((bound_b + bound_a) / 2);
    profit = calc_profit_Sperry_ci_one_line(opt_ci);

  
    return;

  }


double Leaf::optimise_psi_stem_Sperry_Newton() {
  
  // double psi_stem_next = psi_soil_;

  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    psi_stem_next = psi_soil_;
    return(psi_stem_next);
  }
  // optimise for stem water potential
  double diff_value = 0.01; 
  double epsilon = 0.001;
  
  double psi_stem_initial;

  int finished=1;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  psi_stem_next = psi_soil_ + diff_value;

  while(finished == 1){

    psi_stem_initial = psi_stem_next;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;

    y_0 = calc_profit_Sperry(psi_stem_initial);
    y_1 = calc_profit_Sperry(x_1);
    y_2 = calc_profit_Sperry(x_2);


    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    psi_stem_next = psi_stem_initial -  first_dev/sec_dev;

    if(abs(psi_stem_next - psi_stem_initial) < epsilon){
      finished = 0;
    }

  }
    return psi_stem_next;
  }


void Leaf::optimise_psi_stem_Sperry_Newton_recall_one_line(double psi_guess) {

  int count = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    opt_psi_stem = psi_soil_;
    y_0 = calc_profit_Sperry_one_line(opt_psi_stem);
    return;
  }
  // optimise for stem water potential
  double diff_value = 0.001; 
  // epsilon_leaf = 0.0001;
  
  double psi_stem_initial;

  int finished=1;

  opt_psi_stem = psi_guess;

  if (R_IsNA(opt_psi_stem) | (opt_psi_stem > psi_crit)){
// std::cout << "psi_soil" << psi_soil_ <<  std::endl;

  opt_psi_stem = psi_crit - diff_value;

  }

  while(finished == 1){

    count += 1;

    psi_stem_initial = opt_psi_stem;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;
    
    y_2 = calc_profit_Sperry_one_line(x_2);

    y_1 = calc_profit_Sperry_one_line(x_1);

    y_0 = calc_profit_Sperry_one_line(psi_stem_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    opt_psi_stem = psi_stem_initial -  first_dev/sec_dev;
    
    if((opt_psi_stem > psi_crit) | (opt_psi_stem < psi_soil_)){
      
      opt_psi_stem = optimise_psi_stem_Sperry_one_line();
      finished = 0;

    // std::cout << " opt_psi_stem_gss: " << opt_psi_stem <<  std::endl;

      profit = calc_profit_Sperry_one_line(opt_psi_stem);

      return;
    }

    if(abs(opt_psi_stem - psi_stem_initial) < epsilon_leaf){
      finished = 0;
    }

    profit = y_0;
  }
    if(count > 1){
    // std::cout << " " << count;
    }
    return;
  }

void Leaf::optimise_ci_Sperry_Newton_recall_one_line(double ci_guess) {

  int count = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  // Early exit -- XXXX 
  if (psi_soil_ > psi_crit){

    // std::cout << "psi_soil_exceeded" << std::endl;

    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0.0;
    E = 0.0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.001; 
  // TODO: rename as delta_ci
  
  // epsilon_leaf = 0.001;
  
  double ci_initial;

  int unfinished=1; // TODO rename unfinished

  opt_ci = ci_guess;

// add in the max ci
  if (R_IsNA(opt_ci)){
  
  opt_ci = find_max_ci_one_line() - diff_value;

  }

  while(unfinished == 1){

// std::cout << "opt_ci" << opt_ci << "ci_guess" << ci_guess << "ci_initial" << ci_initial << "PPFD" << PPFD_ << "k_l_max" << k_l_max_ << "psi_crit" << psi_crit << "psi_soil" << psi_soil_ << std::endl;


    count += 1;

    if(count > 2){
      
      double max_ci = find_max_ci_one_line();
      optimise_ci_Sperry_one_line(max_ci);

      return;

    }

    ci_initial = opt_ci;

    if (!util::is_finite(opt_ci)) {
      util::stop("Detected NAN ci value");
    }

    x_1 = ci_initial - diff_value;
    x_2 = ci_initial + diff_value;

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream
    double benefit_ = calc_A_lim_one_line(x_2);


    double g_c_ci = (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - x_2); 
    

    E = g_c_ci * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;  
//move to top

    double E_max = calc_E_supply(psi_crit);

  // std::cout << "Benefit " << benefit_ << "E " << E << "E_max " << E_max << "psi_soil" << psi_soil_  << "PPFD "<< PPFD_ << std::endl;

    if(((E_max < E) & (abs(E - E_max) > 1e-15)) | (E < 0)){
    //E_max already known, no need to do find_max_ci_one_line_again 
    opt_ci = find_max_ci_one_line() - diff_value;

    continue;

    } else{

    // std::cout << "E" << E << std::endl;  
    double psi_stem = calc_psi_stem_ci(E);

    double cost_ = calc_hydraulic_cost_Sperry(psi_stem);
    y_2 = benefit_ - lambda_*cost_;
    }

    y_1 = calc_profit_Sperry_ci_one_line(x_1);

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream

    y_0 = calc_profit_Sperry_ci_one_line(ci_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);


//added this in to account for situations where profit curve is flat and 0, likely due to no light
    if(y_0 == 0.0 & first_dev == 0.0 & sec_dev == 0.0){
    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0.0;
    E = 0.0;

    std::cout << "no profit curve" << std::endl;

    return;
    }

    opt_ci = ci_initial -  first_dev/sec_dev;

    if(abs(opt_ci - ci_initial) < (ci_initial*epsilon_leaf)){
     
      unfinished = 0;

    }

    profit = y_0;
  }

    return;
  }

void Leaf::optimise_ci_Sperry_Newton_recall_one_line_max(double ci_guess) {

  int count = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  if (psi_soil_ > psi_crit){
    opt_ci = gamma_25*umol_per_mol_to_Pa;
    profit = 0;
    E = 0;
    return;
  }

  // optimise for stem water potential
  double diff_value = 0.001; 
  // epsilon_leaf = 0.0001;
  
  double ci_initial;

  int finished=1;

  opt_ci = ci_guess;

double max_ci =find_max_ci_one_line();


  while(finished == 1){


// add in the max ci
  if (R_IsNA(opt_ci) | (opt_ci > (max_ci - diff_value))){
  
  opt_ci = max_ci - diff_value;

  }

    count += 1;
    if(count > 0){


      if(count > 15){
      util::stop("stuck");
    }
    }

    ci_initial = opt_ci;

    if (!util::is_finite(opt_ci)) {



      util::stop("Detected NAN ci value");
    }

    x_1 = ci_initial - diff_value;
    x_2 = ci_initial + diff_value;

    y_2 = calc_profit_Sperry_ci_one_line(x_2);


    y_1 = calc_profit_Sperry_ci_one_line(x_1);

// start with calculation of highest ci, this is the estimate with the highest required transpiration stream so the estimate most likely to exceed maximum possible transpiration stream

    y_0 = calc_profit_Sperry_ci_one_line(ci_initial);

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    opt_ci = ci_initial -  first_dev/sec_dev;
      // std::cout <<  "diff" << ci_initial*epsilon;

    if(abs(opt_ci - ci_initial) < (ci_initial*epsilon_leaf)){
      finished = 0;
    }

    profit = y_0;
  }

  // std::cout << "eps" << epsilon_leaf << "diff" << abs(opt_ci - ci_initial) << "tolerated_diff" << ci_initial*epsilon_leaf << "opt_ci" << opt_ci<< "ci_initial" << ci_initial;

    return;
  }

void Leaf::optimise_psi_stem_Sperry_Newton_recall_one_line_pass() {

  int count = 0, errors=0;

  // double psi_stem_next = psi_soil_;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    opt_psi_stem = psi_soil_;
    profit = 0;
    return;
  }
  // optimise for stem water potential
  double diff_value = 0.001; 
  double epsilon = 0.01;
  
  double psi_stem_initial;

  int finished=1;
  int gss = 0;

  // opt_psi_stem = psi_guess;


  if (R_IsNA(opt_psi_stem) | (opt_psi_stem > psi_crit)){

// <<  "psi_crit" << psi_crit 

  opt_psi_stem = psi_soil_ + diff_value;

  // std::cout << " opt_psi_stem_reset: " << opt_psi_stem; 

  }
    // std::cout << "psi_stem_next" << psi_stem_next <<std::endl;

  while(finished == 1){

    count += 1;

    psi_stem_initial = opt_psi_stem;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;
    
    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "psi_stem_initial" << psi_stem_initial << "PPFD_" << PPFD_ << "psi_soil_" << psi_soil_ << "k_l_max_" << k_l_max_ << "psi_crit" << psi_crit << std::endl;

    
    y_0 = calc_profit_Sperry_one_line(psi_stem_initial);
    
    // std::cout << "y_0" <<  y_0 <<std::endl;

    y_1 = calc_profit_Sperry_one_line(x_1);
    
    // std::cout << "y_1" <<  y_1 <<std::endl;

    y_2 = calc_profit_Sperry_one_line(x_2);

    // std::cout << "y_2" <<  y_2 <<std::endl;

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    opt_psi_stem = psi_stem_initial -  first_dev/sec_dev;

    if((opt_psi_stem > psi_crit) | (opt_psi_stem < psi_soil_)){
      
      // std::cout << "errors: " << errors;

      if(errors < 2) {
        opt_psi_stem = std::max(psi_soil_, 0.5*(psi_stem_initial + psi_soil_));
        errors+=1;
      }  else {

        // std::cout << "gss_psi_soil: " << psi_soil_<< " PPFD: " << PPFD_  << " count: " << count << " opt_psi_stem_guess: " << psi_guess << "opt_psi_stem: " << opt_psi_stem;
        // std::cout  <<  " gss: " << psi_stem_initial <<  " opt_psi_stem_updated: "  << opt_psi_stem;

      opt_psi_stem = optimise_psi_stem_Sperry_one_line();
      gss = 1;
      finished = 0;

    // std::cout << " opt_psi_stem_gss: " << opt_psi_stem <<  std::endl;

      profit = calc_profit_Sperry_one_line(opt_psi_stem);

      return;
    }
    }

    if(abs(opt_psi_stem - psi_stem_initial) < epsilon){
      finished = 0;
    }

    profit = y_0;
  }
    if(count > 1){
    // std::cout << " " << count;
    }
    return;
  }

double Leaf::optimise_psi_stem_Sperry_Newton_recall() {
  
  // double psi_stem_next = psi_soil_;

  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    psi_stem_next = psi_soil_;
    return(psi_stem_next);
  }
  // optimise for stem water potential
  double diff_value = 0.001; 
  double epsilon = 0.001;
  
  double psi_stem_initial;

  int finished=1;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

    // std::cout << "psi_stem_next" << psi_stem_next <<std::endl;

  if (R_IsNA(psi_stem_next)){
  psi_stem_next = psi_soil_ + (psi_crit - psi_soil_)/2;
  }

      // std::cout << "psi_stem_next" << psi_stem_next <<std::endl;

  while(finished == 1){

    psi_stem_initial = psi_stem_next;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;

    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "psi_stem_initial" << psi_stem_initial <<std::endl;


    y_0 = calc_profit_Sperry(psi_stem_initial);

    // std::cout << "y_0" <<  y_0 <<std::endl;

    y_1 = calc_profit_Sperry(x_1);

    // std::cout << "y_1" <<  y_1 <<std::endl;

    y_2 = calc_profit_Sperry(x_2);
    
    // std::cout << "y_2" <<  y_2 <<std::endl;


    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    psi_stem_next = psi_stem_initial -  first_dev/sec_dev;

    if(abs(psi_stem_next - psi_stem_initial) < epsilon){
      finished = 0;
    }

  }
    return psi_stem_next;
  }

double Leaf::optimise_ci_Sperry_Newton() {

// double psi_stem_next = psi_soil_;

  if ((PPFD_ < 1.5e-8 )| (psi_soil_ > psi_crit)){
    c_i_next = gamma_25*umol_per_mol_to_Pa;
    return(c_i_next);
  }
  // optimise for stem water potential
  double diff_value = 0.1; 
  double epsilon = 0.001;
  
  double c_i_initial;

  int finished=1;

  // double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  c_i_next = ((gamma_25*umol_per_mol_to_Pa + diff_value) + find_max_ci())/2;

    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;

  while(finished == 1){

    c_i_initial = c_i_next;

 double x_1 = c_i_initial - diff_value;
 double x_2 = c_i_initial + diff_value;

    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;


  double y_0 = calc_profit_Sperry_ci(c_i_initial);
  
    // std::cout << "y_0" <<  y_0 <<std::endl;

  double y_1 = calc_profit_Sperry_ci(x_1);
  
    // std::cout <<"y_1" << y_1  << std::endl;

  double y_2 = calc_profit_Sperry_ci(x_2);

  // std::cout  << "y_2" << y_2 << std::endl;

  double first_dev = (y_2 - y_1)/(2*diff_value);
  double sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    c_i_next = c_i_initial -  first_dev/sec_dev;

    if(abs(c_i_next - c_i_initial) < epsilon){
      finished = 0;
    }

  }
    return c_i_next;
  }


double Leaf::optimise_ci_Sperry_Newton_recall() {

// double psi_stem_next = psi_soil_;

  if ((PPFD_ < 1.5e-8) | (psi_soil_ > psi_crit)){
    c_i_next = gamma_25*umol_per_mol_to_Pa;
    return(c_i_next);
  }
  // optimise for stem water potential
  double diff_value = 0.01; 
  double epsilon = 0.001;
  
  double c_i_initial;

  int finished=1;

  // double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  if (R_IsNA(c_i_next)){
  c_i_next = ((gamma_25*umol_per_mol_to_Pa + diff_value) + find_max_ci())/2;
  }

  
    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;

  while(finished == 1){

    c_i_initial = c_i_next;

 double x_1 = c_i_initial - diff_value;
 double x_2 = c_i_initial + diff_value;

    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;


  double y_0 = calc_profit_Sperry_ci(c_i_initial);
  
    // std::cout << "y_0" <<  y_0 <<std::endl;

  double y_1 = calc_profit_Sperry_ci(x_1);
  
    // std::cout <<"y_1" << y_1  << std::endl;

  double y_2 = calc_profit_Sperry_ci(x_2);

  // std::cout  << "y_2" << y_2 << std::endl;

  double first_dev = (y_2 - y_1)/(2*diff_value);
  double sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    c_i_next = c_i_initial -  first_dev/sec_dev;

    if(abs(c_i_next - c_i_initial) < epsilon){
      finished = 0;
    }

  }
    return c_i_next;
  }

// Bartlett et al. implementation. Predicting shifts in the functional composition of tropical forests under increased drought and CO2 from trade-offs among plant hydraulic traits. Bartlett et al. (2018).

double Leaf::calc_hydraulic_cost_Bartlett(double psi_stem) {

  double k_l_soil_ = k_l_max_ * calc_cond_vuln(psi_soil_);
  double k_l_stem_ = k_l_max_ * calc_cond_vuln(psi_stem);
  double height_ = (K_s * huber_value)/k_l_max_;
  // std::cout << "K_s_ " << K_s << "huber_value " << huber_value << std::endl;

  return beta * huber_value * height_ * pow((1 - k_l_stem_ / k_l_soil_), beta_2);
}

double Leaf::calc_profit_Bartlett(double psi_stem) {                                  

  double benefit_ =
      calc_assim_gross(psi_stem);

  double cost_ =
      calc_hydraulic_cost_Bartlett(psi_stem);

  return benefit_ - cost_;
}

// need docs on Golden Section Search and reference to Bartlett.
double Leaf::optimise_psi_stem_Bartlett() {

  double gr = (sqrt(5) + 1) / 2;
  double opt_psi_stem = psi_soil_;

  // optimise for stem water potential
  if (psi_soil_ < psi_crit) {
    double bound_a = psi_soil_;
    double bound_b = psi_crit;

    double delta_crit = 1e-3;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Bartlett(bound_c);

      double profit_at_d =
          calc_profit_Bartlett(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem = ((bound_b + bound_a) / 2);
  }
  // std::cout << "psi_crit " << psi_crit << std::endl;

  return opt_psi_stem;
}



double Leaf::optimise_ci_Bartlett() {

 double gr = (sqrt(5) + 1) / 2;

  // optimise for stem water potential
    double bound_a = 0;
    double bound_b = ca;

    double delta_crit = 1e-3;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Bartlett_ci(bound_c);

      double profit_at_d =
          calc_profit_Bartlett_ci(bound_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

   double opt_ci = ((bound_b + bound_a) / 2);

  // std::cout << "psi_crit " << psi_crit << std::endl;

  return opt_ci;
}

double Leaf::calc_profit_Bartlett_ci(double c_i) {                                  

  double benefit_ =
      calc_assim_gross_ci(c_i);
  
  double g_c = (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - c_i); 
  // std::cout << "g_c" << g_c << std::endl;
  double E_ci = g_c * 1.6 * atm_vpd_ / kg_to_mol_h2o / atm_kpa;
  double psi_stem = calc_psi_stem_ci(E_ci);
  double cost_ = calc_hydraulic_cost_Bartlett(psi_stem);

  // std::cout << "kPa_to_Pa"<<kPa_to_Pa << "atm_kpa" << atm_kpa << "umol_per_mol_to_mol_per_mol" <<umol_per_mol_to_mol_per_mol<< "benefit" << benefit_ << "g_c" << g_c <<"E" << E_ci <<"psi_stem" << psi_stem << "cost_" << cost_  << "ci" << c_i <<std::endl;


  return benefit_ - cost_;
}

double Leaf::calc_assim_gross_ci(double ci) {

  A_lim = calc_A_lim(ci);
  return A_lim;
  
}

} // namespace plant