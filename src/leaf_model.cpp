#include <plant/leaf_model.h>
#include <cmath>

namespace plant {
Leaf::Leaf(double vcmax, double p_50, double c, double b,
           double psi_crit, // derived from b and c
           double beta, double beta_2, double huber_value, double K_s)
    : vcmax(vcmax), p_50(p_50), c(c), b(b), psi_crit(psi_crit), beta(beta),
      beta_2(beta_2), huber_value(huber_value), K_s(K_s), ci(NA_REAL), g_c(NA_REAL),
      A_lim(NA_REAL), E(NA_REAL), psi(NA_REAL), profit(NA_REAL), psi_stem_next(NA_REAL), c_i_next(NA_REAL), lambda_(NA_REAL),
      PPFD_(NA_REAL), psi_soil_(NA_REAL), k_l_max_(NA_REAL) {
    setup_E_supply(100);
    // setup_psi(100);
}

void Leaf::set_physiology(double PPFD, double psi_soil, double k_l_max) {
   lambda_ = calc_assim_gross(PPFD, psi_soil, psi_crit, k_l_max) / (k_l_max * calc_cond_vuln(psi_soil) - k_l_max * calc_cond_vuln(psi_crit));
   PPFD_ = PPFD;
   psi_soil_ = psi_soil;
   k_l_max_ = k_l_max;
}

double Leaf::p_50_to_K_s() const {
  return std::pow(1.638999 * (p_50 / 2), -1.38);
}

double Leaf::calc_vul_b() const {
  return p_50 / pow(-log(1 - 50.0 / 100.0), 1 / c);
}

// REMOVED k_l_max
// integrates, returns conductivity at given psi_stem kg m^-2 s^-1 MPa^-1

double Leaf::calc_cond_vuln(double psi) const {
  return exp(-pow((psi / b), c));
}

// REMOVED k_l_max
void Leaf::setup_E_supply(double resolution) {
    // integrate and accumulate results
    auto x_psi = std::vector<double>{0.0};  // {0.0}
    auto y_cumulative_E = std::vector<double>{0.0}; // {0.0}
    double step = (psi_crit+1)/resolution;
   
    for (double psi_spline = 0 + step; psi_spline <= (psi_crit + 1); psi_spline += step) {
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
   double Leaf::calc_E_supply_full_integration(double k_l_max, double psi_soil,
                              double psi_stem) {
 std::function<double(double)> f;
 f = [&](double psi) -> double { return calc_cond_vuln(psi); };

 return k_l_max * integrator.integrate(f, psi_soil, psi_stem);
                              }

double Leaf::calc_E_supply(double psi_stem) {
    // integration of calc_cond_vuln over [psi_soil, psi_stem]
    return k_l_max_ * (E_from_psi.eval(psi_stem) - E_from_psi.eval(psi_soil_));
}

double Leaf::calc_psi_stem_ci(double k_l_max, double psi_soil,
                           double E_ci) {
    // integration of calc_cond_vuln over [psi_soil, psi_stem]
    double E_psi_stem = E_ci/k_l_max +  E_from_psi.eval(psi_soil);
    return psi_from_E.eval(E_psi_stem);
    }

double Leaf::calc_psi_from_E(double E_psi_stem) {
    // integration of calc_cond_vuln over [psi_soil, psi_stem]
    return psi_from_E.eval(E_psi_stem);
    }


double Leaf::calc_min_psi(double PPFD, double psi_soil, double k_l_max) {                                  

  // not clear what x is here
  auto target = [&](double x_ci) mutable -> double {
    return min_psi(PPFD, x_ci, psi_soil, k_l_max);
  };

ci = util::uniroot(target, 0, ca, 1e-6, 1000);
return ci;
}

double Leaf::min_psi(double PPFD, double x_ci, double psi_soil, double k_l_max) {                                  
  double benefit_ =
      calc_A_lim(PPFD, x_ci);

  return (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - x_ci);

}

// returns E kg m^-2 s^-1
double Leaf::calc_g_c(double psi_soil, double psi_stem, double k_l_max) {
  double E_supply = calc_E_supply(psi_stem);

  return atm_kpa * E_supply * kg_to_mol_h2o / atm_vpd / 1.6;
}


double Leaf::calc_j(double PPFD) {
  double jmax = vcmax * vcmax_25_to_jmax_25;
  double j = (a * PPFD + jmax -
              sqrt(pow(a * PPFD + jmax, 2) - 4 * curv_fact * a * PPFD * jmax)) /
             (2 * curv_fact); // check brackets are correct

  return j;           

}

// returns A_c umol m^-2 s^-1
double Leaf::calc_A_c(double ci_) {
  return (vcmax * (ci_ - gamma_25 * umol_per_mol_to_Pa)) / (ci_ + km_25);
}

// returns A_j umol m^-2 s^-1
double Leaf::calc_A_j(double PPFD, double ci_) {

  double j = calc_j(PPFD);
  return j / 4 *
         ((ci_ - gamma_25 * umol_per_mol_to_Pa) /
          (ci_ + 2 * gamma_25 * umol_per_mol_to_Pa));
}

// returns co-limited assimilation umol m^-2 s^-1
double Leaf::calc_A_lim(double PPFD, double ci_) {

  double A_c = calc_A_c(ci_);
  double A_j = calc_A_j(PPFD, ci_);
  // double R_d = vcmax * 0.015;

  return (A_c + A_j - sqrt(pow(A_c + A_j, 2) - 4 * 0.98 * A_c * A_j)) /
             (2 * 0.98);
}


// returns co-limited assimilation umol m^-2 s^-1
double Leaf::calc_A_lim_one_line(double PPFD, double ci_) {

  double c2 = 13.13652;

  double j = calc_j(PPFD);
  return j / 4 *
         ((ci_ - gamma_25 * umol_per_mol_to_Pa) /
          (ci_ + c2));
}



// returns difference between co-limited assimilation and g_c, to be minimised (umol m^-2 s^-1)
double Leaf::diff_ci(double PPFD, double x, double psi_soil, double psi_stem, double k_l_max) {

  double A_lim_ = calc_A_lim(PPFD, x);

  double g_c_ = calc_g_c(psi_soil, psi_stem, k_l_max);

  return A_lim_ * umol_per_mol_to_mol_per_mol -
         (g_c_ * (ca - x) / (atm_kpa * kPa_to_Pa));
}

// need to fill in tol and max_iteratiosn
double Leaf::calc_assim_gross(
    double PPFD, double psi_soil, double psi_stem, double k_l_max) {

 if (psi_soil == psi_stem){
    ci = gamma_25*umol_per_mol_to_Pa;
      } else{


  // not clear what x is here
  auto target = [&](double x) mutable -> double {
    return diff_ci(PPFD, x, psi_soil, psi_stem, k_l_max);
  };

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  ci = util::uniroot(target, 0, ca, 1e-6, 1000);
      }
  A_lim = calc_A_lim(PPFD, ci);

  g_c = calc_g_c(psi_soil, psi_stem, k_l_max);

  //E is in m^3 m-2 

  E = g_c * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;

  psi = psi_stem;

  return A_lim;
  
}

// Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::calc_hydraulic_cost_Sperry(double psi_soil, double psi_stem,
                                 double k_l_max) {
  double k_l_soil_ = k_l_max * calc_cond_vuln(psi_soil);
  double k_l_stem_ = k_l_max * calc_cond_vuln(psi_stem);

  return k_l_soil_ - k_l_stem_;
}

double Leaf::calc_lambda_Sperry(double PPFD, double psi_soil, double k_l_max){
  double max_assim = calc_assim_gross(PPFD, psi_soil, psi_crit, k_l_max);
  double max_cost = k_l_max * calc_cond_vuln(psi_soil) - k_l_max * calc_cond_vuln(psi_crit);
  return max_assim/max_cost;
}

double Leaf::find_max_ci(
    double PPFD, double psi_soil, double k_l_max) {


  // not clear what x is here
  auto target = [&](double x) mutable -> double {
    return diff_ci(PPFD, x, psi_soil, psi_crit, k_l_max);
  };

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  ci = util::uniroot(target, 0, ca, 1e-6, 1000);

  return ci;
  
}


double Leaf::calc_profit_Sperry(double PPFD, double psi_soil, double psi_stem, double k_l_max) {

  double benefit_ =
      calc_assim_gross(PPFD, psi_soil, psi_stem, k_l_max);
  double cost_ = calc_hydraulic_cost_Sperry(psi_soil, psi_stem, k_l_max);

  return benefit_ - lambda_ * cost_;
}





// need docs on Golden Section Search and reference to Bartlett.
double Leaf::optimise_psi_stem_Sperry(double PPFD, double psi_soil, double k_l_max) {

  double gr = (sqrt(5) + 1) / 2;
  double opt_psi_stem = psi_soil;


  if ((PPFD < 1.5e-8 )| (psi_soil > psi_crit)){
    return(opt_psi_stem);
  }

  // optimise for stem water potential
    double bound_a = psi_soil;
    double bound_b = psi_crit;

    double delta_crit = 1e-3;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Sperry(PPFD, psi_soil, bound_c, k_l_max);

      double profit_at_d =
          calc_profit_Sperry(PPFD, psi_soil, bound_d, k_l_max);

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
  

double Leaf::optimise_psi_stem_Sperry_one_line(double PPFD, double psi_soil, double k_l_max) {

  double gr = (sqrt(5) + 1) / 2;
  double opt_psi_stem = psi_soil;


  if ((PPFD < 1.5e-8 )| (psi_soil > psi_crit)){
    return(opt_psi_stem);
  }

  // optimise for stem water potential
    double bound_a = psi_soil;
    double bound_b = psi_crit;

    double delta_crit = 1e-3;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Sperry_one_line(PPFD, psi_soil, bound_c, k_l_max);

      double profit_at_d =
          calc_profit_Sperry_one_line(PPFD, psi_soil, bound_d, k_l_max);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem = ((bound_b + bound_a) / 2);

    std::cout << "hello" << std::endl;
  


    return opt_psi_stem;

  }


double Leaf::optimise_psi_stem_Sperry_Newton(double PPFD, double psi_soil, double k_l_max) {
  
  // double psi_stem_next = psi_soil;

  if ((PPFD < 1.5e-8 )| (psi_soil > psi_crit)){
    psi_stem_next = psi_soil;
    return(psi_stem_next);
  }
  // optimise for stem water potential
  double diff_value = 0.01; 
  double epsilon = 0.001;
  
  double psi_stem_initial;

  int finished=1;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  psi_stem_next = psi_soil + diff_value;

  while(finished == 1){

    psi_stem_initial = psi_stem_next;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;

    y_0 = calc_profit_Sperry(PPFD, psi_soil, psi_stem_initial, k_l_max);
    y_1 = calc_profit_Sperry(PPFD, psi_soil, x_1, k_l_max);
    y_2 = calc_profit_Sperry(PPFD, psi_soil, x_2, k_l_max);


    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    psi_stem_next = psi_stem_initial -  first_dev/sec_dev;

    if(abs(psi_stem_next - psi_stem_initial) < epsilon){
      finished = 0;
    }

  }
    return psi_stem_next;
  }

double Leaf::calc_assim_gross_one_line(double PPFD, double psi_soil, double psi_stem, double k_l_max) {

    if (psi_soil == psi_stem){
    g_c = 0;
    ci = gamma_25*umol_per_mol_to_Pa;
      } else{
    g_c = calc_g_c(psi_soil, psi_stem, k_l_max);
    
    double c2 = 13.13652;
    
    double first_term = 8 * calc_j(PPFD = PPFD) * umol_per_mol_to_mol_per_mol * (atm_kpa*kPa_to_Pa) * g_c * (-ca + c2 + 2 * gamma_25 * umol_per_mol_to_Pa);
    double second_term = 16 * pow(g_c, 2);
    double third_term = pow((ca + c2),2);
    double fourth_term = pow(calc_j(PPFD = PPFD), 2) * pow(umol_per_mol_to_mol_per_mol,2) * pow(atm_kpa*kPa_to_Pa, 2);
    double fifth_term = 4*ca*g_c;
    double sixth_term = 4*c2*g_c;
    double seventh_term = calc_j(PPFD = PPFD)*umol_per_mol_to_mol_per_mol*(atm_kpa*kPa_to_Pa);
    double eigth_term = 8*g_c;

    ci = (sqrt(first_term + second_term*third_term+ fourth_term) + fifth_term - sixth_term - seventh_term)/eigth_term;
      } 
    A_lim = calc_A_lim_one_line(PPFD, ci);
    
    E = g_c * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;
    
    psi = psi_stem;
    
    return A_lim;
}

double Leaf::calc_profit_Sperry_one_line(double PPFD, double psi_soil, double psi_stem, double k_l_max) {

  double benefit_ =
      calc_assim_gross_one_line(PPFD, psi_soil, psi_stem, k_l_max);
  double cost_ = calc_hydraulic_cost_Sperry(psi_soil, psi_stem, k_l_max);

  return benefit_ - lambda_ * cost_;
}


double Leaf::optimise_psi_stem_Sperry_Newton_recall_one_line(double PPFD, double psi_soil, double k_l_max) {
  
  // double psi_stem_next = psi_soil;

  if ((PPFD < 1.5e-8 )| (psi_soil > psi_crit)){
    psi_stem_next = psi_soil;
    return(psi_stem_next);
  }
  // optimise for stem water potential
  double diff_value = 0.01; 
  double epsilon = 0.001;
  
  double psi_stem_initial;

  int finished=1;
  int gss = 0;

  double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

    // std::cout << "psi_stem_next" << psi_stem_next <<std::endl;


  if (R_IsNA(psi_stem_next)){
  psi_stem_next = psi_soil + diff_value;
  }
    // std::cout << "psi_stem_next" << psi_stem_next <<std::endl;

  while(finished == 1){

    psi_stem_initial = psi_stem_next;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;
    
    std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "psi_stem_initial" << psi_stem_initial << "PPFD" << PPFD << "psi_soil" << psi_soil << "k_l_max" << k_l_max << "psi_crit" << psi_crit << std::endl;

    
    y_0 = calc_profit_Sperry_one_line(PPFD, psi_soil, psi_stem_initial, k_l_max);
    
    // std::cout << "y_0" <<  y_0 <<std::endl;

    y_1 = calc_profit_Sperry_one_line(PPFD, psi_soil, x_1, k_l_max);
    
    // std::cout << "y_1" <<  y_1 <<std::endl;

    y_2 = calc_profit_Sperry_one_line(PPFD, psi_soil, x_2, k_l_max);

    // std::cout << "y_2" <<  y_2 <<std::endl;

    first_dev = (y_2 - y_1)/(2*diff_value);
    sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    psi_stem_next = psi_stem_initial -  first_dev/sec_dev;

    if(psi_stem_next > psi_crit | psi_stem_next < psi_soil){
      psi_stem_next = optimise_psi_stem_Sperry_one_line(PPFD, psi_soil, k_l_max);
      gss = 1;
      finished = 0;
    }

    if(abs(psi_stem_next - psi_stem_initial) < epsilon){
      finished = 0;
    }

  }
    std::cout << "psi_stem_next" << psi_stem_next << "gss" << gss <<  std::endl;

    return psi_stem_next;
  }

double Leaf::optimise_psi_stem_Sperry_Newton_recall(double PPFD, double psi_soil, double k_l_max) {
  
  // double psi_stem_next = psi_soil;

  if ((PPFD < 1.5e-8 )| (psi_soil > psi_crit)){
    psi_stem_next = psi_soil;
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
  psi_stem_next = psi_soil + (psi_crit - psi_soil)/2;
  }

      // std::cout << "psi_stem_next" << psi_stem_next <<std::endl;

  while(finished == 1){

    psi_stem_initial = psi_stem_next;

    x_1 = psi_stem_initial - diff_value;
    x_2 = psi_stem_initial + diff_value;

    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "psi_stem_initial" << psi_stem_initial <<std::endl;


    y_0 = calc_profit_Sperry(PPFD, psi_soil, psi_stem_initial, k_l_max);

    // std::cout << "y_0" <<  y_0 <<std::endl;

    y_1 = calc_profit_Sperry(PPFD, psi_soil, x_1, k_l_max);

    // std::cout << "y_1" <<  y_1 <<std::endl;

    y_2 = calc_profit_Sperry(PPFD, psi_soil, x_2, k_l_max);
    
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

  







double Leaf::optimise_ci_Sperry_Newton(double PPFD, double psi_soil, double k_l_max) {

// double psi_stem_next = psi_soil;

  if ((PPFD < 1.5e-8 )| (psi_soil > psi_crit)){
    c_i_next = gamma_25*umol_per_mol_to_Pa;
    return(c_i_next);
  }
  // optimise for stem water potential
  double diff_value = 0.1; 
  double epsilon = 0.001;
  
  double c_i_initial;

  int finished=1;

  // double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;

  c_i_next = ((gamma_25*umol_per_mol_to_Pa + diff_value) + find_max_ci(PPFD, psi_soil, k_l_max))/2;

    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;

  while(finished == 1){

    c_i_initial = c_i_next;

 double x_1 = c_i_initial - diff_value;
 double x_2 = c_i_initial + diff_value;

    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;


  double y_0 = calc_profit_Sperry_ci(PPFD, psi_soil, c_i_initial, k_l_max);
  
    // std::cout << "y_0" <<  y_0 <<std::endl;

  double y_1 = calc_profit_Sperry_ci(PPFD, psi_soil, x_1, k_l_max);
  
    // std::cout <<"y_1" << y_1  << std::endl;

  double y_2 = calc_profit_Sperry_ci(PPFD, psi_soil, x_2, k_l_max);

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


double Leaf::optimise_ci_Sperry_Newton_recall(double PPFD, double psi_soil, double k_l_max) {

// double psi_stem_next = psi_soil;

  if (PPFD < 1.5e-8 | psi_soil > psi_crit){
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
  c_i_next = ((gamma_25*umol_per_mol_to_Pa + diff_value) + find_max_ci(PPFD, psi_soil, k_l_max))/2;
  }

  
    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;

  while(finished == 1){

    c_i_initial = c_i_next;

 double x_1 = c_i_initial - diff_value;
 double x_2 = c_i_initial + diff_value;

    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;


  double y_0 = calc_profit_Sperry_ci(PPFD, psi_soil, c_i_initial, k_l_max);
  
    // std::cout << "y_0" <<  y_0 <<std::endl;

  double y_1 = calc_profit_Sperry_ci(PPFD, psi_soil, x_1, k_l_max);
  
    // std::cout <<"y_1" << y_1  << std::endl;

  double y_2 = calc_profit_Sperry_ci(PPFD, psi_soil, x_2, k_l_max);

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



double Leaf::optimise_ci_Sperry_Newton_recall_one_line(double PPFD, double psi_soil, double k_l_max) {

// double psi_stem_next = psi_soil;

  if (PPFD < 1.5e-8 | psi_soil > psi_crit){
    c_i_next = gamma_25*umol_per_mol_to_Pa;
    return(c_i_next);
  }
  // optimise for stem water potential
  double diff_value = 0.01; 
  double epsilon = 0.001;
  
  double c_i_initial;

  int finished=1;

  // double x_1, x_2, y_0, y_1, y_2, first_dev, sec_dev;
  
  double max_lim =  calc_assim_gross_one_line(PPFD, psi_soil, psi_crit, k_l_max);

  double max_ci = ci;

  std::cout << "max_ci" << max_ci <<std::endl;


  if (R_IsNA(c_i_next)){
  c_i_next = ((gamma_25*umol_per_mol_to_Pa + diff_value) + max_ci)/2;
  }
  
  if (c_i_next > max_ci - diff_value){
  c_i_next = ((gamma_25*umol_per_mol_to_Pa + diff_value) + max_ci)/2;
  }
    // std::cout << "x_1" << x_1 << "x_2" << x_2 << "diff_value" << diff_value << "c_i_initial" << c_i_initial <<std::endl;

  while(finished == 1){

    c_i_initial = c_i_next;

 double x_1 = c_i_initial - diff_value;
 double x_2 = c_i_initial + diff_value;

    std::cout << "c_i_initial" << c_i_initial << "c_i_next" << c_i_next << "x_1" << x_1 << "x_2" << x_2 << "psi_soil" << psi_soil << "k_l_max" << k_l_max << "PPFD" << PPFD <<std::endl;


  double y_0 = calc_profit_Sperry_ci_one_line(PPFD, psi_soil, c_i_initial, k_l_max);
  
    // std::cout << "y_0" <<  y_0 <<std::endl;

  double y_1 = calc_profit_Sperry_ci_one_line(PPFD, psi_soil, x_1, k_l_max);
  
    // std::cout <<"y_1" << y_1  << std::endl;

  double y_2 = calc_profit_Sperry_ci_one_line(PPFD, psi_soil, x_2, k_l_max);

  // std::cout  << "y_2" << y_2 << std::endl;

  double first_dev = (y_2 - y_1)/(2*diff_value);
  double sec_dev = (y_2 - 2*y_0 + y_1)/pow(diff_value, 2);

    c_i_next = c_i_initial -  first_dev/sec_dev;

    if(abs(c_i_next - c_i_initial) < epsilon){
      finished = 0;
    }

  }

    std::cout << "c_i_next" << c_i_next <<std::endl;

    return c_i_next;
  }

double Leaf::calc_profit_Sperry_ci(double PPFD, double psi_soil, double c_i,double k_l_max) {                                  
  double benefit_ =
      calc_A_lim(PPFD, c_i);
  double g_c_ci = (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - c_i); 
  double E_ci = g_c_ci * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;
  // std::cout << "E_ci" << E_ci << "g_c" << g_c_ci << "c_i" << c_i << std::endl;
  
  double psi_stem = calc_psi_stem_ci(k_l_max, psi_soil, E_ci);
  double cost_ = calc_hydraulic_cost_Sperry(psi_soil, psi_stem, k_l_max);
  double lambda_ = calc_assim_gross(PPFD, psi_soil, psi_crit, k_l_max) / calc_hydraulic_cost_Sperry(psi_soil, psi_crit, k_l_max);

  return benefit_ - lambda_*cost_;
}


double Leaf::calc_profit_Sperry_ci_one_line(double PPFD, double psi_soil, double c_i,double k_l_max) {                                  
  double benefit_ =
      calc_A_lim_one_line(PPFD, c_i);
  double g_c_ci = (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - c_i); 
  double E_ci = g_c_ci * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;
  // std::cout << "E_ci" << E_ci << "g_c" << g_c_ci << "c_i" << c_i << std::endl;
  
  double psi_stem = calc_psi_stem_ci(k_l_max, psi_soil, E_ci);
  double cost_ = calc_hydraulic_cost_Sperry(psi_soil, psi_stem, k_l_max);

  return benefit_ - lambda_*cost_;
}


// double Leaf::calc_profit_Sperry_ci(double PPFD, double psi_soil, double c_i,double k_l_max) {                                  
  
//   double benefit_ =
//       calc_A_lim(PPFD, c_i);
//   double g_c_ci = (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - c_i); 
//   double E_ci = g_c_ci * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;
//   std::cout << "E_ci" << E_ci << "g_c" << g_c_ci << "c_i" << c_i << std::endl;
  
//   double psi_stem = calc_psi_stem_ci(k_l_max, psi_soil, E_ci);
//   double cost_ = calc_hydraulic_cost_Sperry(psi_soil, psi_stem, k_l_max);
//   double lambda_ = calc_assim_gross(PPFD, psi_soil, psi_crit, k_l_max) / calc_hydraulic_cost_Sperry(psi_soil, psi_crit, k_l_max);

//   return benefit_ - lambda_*cost_;
// }
// }



































// Bartlett et al. implementation. Predicting shifts in the functional composition of tropical forests under increased drought and CO2 from trade-offs among plant hydraulic traits. Bartlett et al. (2018).

double Leaf::calc_hydraulic_cost_Bartlett(double psi_soil, double psi_stem,
                                          double k_l_max) {

  double k_l_soil_ = k_l_max * calc_cond_vuln(psi_soil);
  double k_l_stem_ = k_l_max * calc_cond_vuln(psi_stem);
  double height_ = (K_s * huber_value)/k_l_max;
  // std::cout << "K_s_ " << K_s << "huber_value " << huber_value << std::endl;

  return beta * huber_value * height_ * pow((1 - k_l_stem_ / k_l_max), beta_2);
}

double Leaf::calc_profit_Bartlett(double PPFD, double psi_soil, double psi_stem,
                                  double k_l_max) {                                  

  double benefit_ =
      calc_assim_gross(PPFD, psi_soil, psi_stem, k_l_max);

  double cost_ =
      calc_hydraulic_cost_Bartlett(psi_soil, psi_stem, k_l_max);

  return benefit_ - cost_;
}

// need docs on Golden Section Search and reference to Bartlett.
double Leaf::optimise_psi_stem_Bartlett(double PPFD, double psi_soil, double k_l_max) {

  double gr = (sqrt(5) + 1) / 2;
  double opt_psi_stem = psi_soil;

  // optimise for stem water potential
  if (psi_soil < psi_crit) {
    double bound_a = psi_soil;
    double bound_b = psi_crit;

    double delta_crit = 1e-3;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Bartlett(PPFD, psi_soil, bound_c, k_l_max);

      double profit_at_d =
          calc_profit_Bartlett(PPFD, psi_soil, bound_d, k_l_max);

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



double Leaf::optimise_ci_Bartlett(double PPFD, double psi_soil, double k_l_max) {

 double gr = (sqrt(5) + 1) / 2;

  // optimise for stem water potential
    double bound_a = 0;
    double bound_b = ca;

    double delta_crit = 1e-3;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Bartlett_ci(PPFD, psi_soil, bound_c, k_l_max);

      double profit_at_d =
          calc_profit_Bartlett_ci(PPFD, psi_soil, bound_d, k_l_max);

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

double Leaf::calc_profit_Bartlett_ci(double PPFD, double psi_soil, double c_i,double k_l_max) {                                  

  double benefit_ =
      calc_assim_gross_ci(PPFD, c_i);
  
  double g_c = (benefit_ * umol_per_mol_to_mol_per_mol * atm_kpa * kPa_to_Pa)/(ca - c_i); 
  // std::cout << "g_c" << g_c << std::endl;
  double E_ci = g_c * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;
  double psi_stem = calc_psi_stem_ci(k_l_max, psi_soil, E_ci);
  double cost_ = calc_hydraulic_cost_Bartlett(psi_soil, psi_stem, k_l_max);

  // std::cout << "kPa_to_Pa"<<kPa_to_Pa << "atm_kpa" << atm_kpa << "umol_per_mol_to_mol_per_mol" <<umol_per_mol_to_mol_per_mol<< "benefit" << benefit_ << "g_c" << g_c <<"E" << E_ci <<"psi_stem" << psi_stem << "cost_" << cost_  << "ci" << c_i <<std::endl;


  return benefit_ - cost_;
}

double Leaf::calc_assim_gross_ci(double PPFD, double ci) {

  A_lim = calc_A_lim(PPFD, ci);
  return A_lim;
  
}

} // namespace plant
