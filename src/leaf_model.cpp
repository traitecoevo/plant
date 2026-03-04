#include <plant/leaf_model.h>
#include <cmath>
#include <plant/models/tf24_environment.h>

namespace plant {
Leaf::Leaf()
    :
    vcmax_25(96), // umol m^-2 s^-1 
    c(2.680147), //unitless
    b(3.898245), //-MPa
    psi_crit(5.870283), //-MPa 
    beta2(1.5), //exponent for effect of hydraulic risk (unitless)
    jmax_25(157.44), // maximum electron transport rate umol m^-2 s^-1
    hk_s(4),  // maximum hydraulic-dependent sapwood turnover rate yr ^ -1
    a(0.30), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(0.7), //curvature factor for the light response curve (unitless)
    curv_fact_colim(0.99), //curvature factor for the colimited photosythnthesis equatiom
    GSS_tol_abs(1e-3),
    vulnerability_curve_ncontrol(100),
    ci_abs_tol(1e-3),
    ci_niter(1000),
    g1_TF24(46.32995) //cost parameter for TF24 profit model umol m^-2 s^-1
   {
      setup_transpiration(100); // arg: num control points for integration
      setup_clean_leaf();
}

Leaf::Leaf(double vcmax_25, double c, double b,
           double psi_crit, // derived from b and c,
           double beta2, double jmax_25, double hk_s,
           double a, double curv_fact_elec_trans, double curv_fact_colim, 
           double GSS_tol_abs,
           double vulnerability_curve_ncontrol,
           double ci_abs_tol,
           double ci_niter,
          double g1_TF24)
    : vcmax_25(vcmax_25), // umol m^-2 s^-1 
    c(c), //unitless
    b(b), //-MPa
    psi_crit(psi_crit), //-MPa 
    beta2(beta2), //exponent for effect of hydraulic risk (unitless)
    jmax_25(jmax_25), // maximum electron transport rate umol m^-2 s^-1
    hk_s(hk_s),  // maximum hydraulic-dependent sapwood turnover rate yr ^ -1
    a(a), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(curv_fact_elec_trans), //curvature factor for the light response curve (unitless)
    curv_fact_colim(curv_fact_colim), //curvature factor for the colimited photosythnthesis equation
    GSS_tol_abs(GSS_tol_abs),
    vulnerability_curve_ncontrol(vulnerability_curve_ncontrol),
    ci_abs_tol(ci_abs_tol),
    ci_niter(ci_niter),
    g1_TF24(g1_TF24) //cost parameter for TF24 profit model umol m^-2 s^-1
   {
      setup_transpiration(vulnerability_curve_ncontrol); // arg: num control points for integration
      setup_clean_leaf();
}

// set various states and physiology parameters obtained from TF24 to NA to clean leaf object
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
  root_collar_psi_ = NA_REAL; //-MPa
  leaf_temp_= NA_REAL; // deg C
  PPFD_= NA_REAL; //umol m^-2 s^-1
  atm_vpd_= NA_REAL; //kPa 
  atm_o2_kpa_= NA_REAL; // kPa
  atm_kpa_= NA_REAL; // kPa
  ca_= NA_REAL; //Pa
  opt_psi_stem_= NA_REAL; //-MPa 
  opt_ci_= NA_REAL; //Pa 
  E_up_ = NA_REAL;
  soil_consumption_.reserve(2);
}

//sets various parameters which are constant for a given node at a given time

void Leaf::set_physiology(double rho, double a_bio, double PPFD, std::vector<double> soil_moist, double root_collar_psi, double leaf_specific_conductance_max, double atm_vpd, double ca, double sapwood_volume_per_leaf_area, double leaf_temp, double atm_o2_kpa, double atm_kpa) {
   rho_ = rho;
   a_bio_ = a_bio;
   atm_vpd_ = atm_vpd;
   leaf_temp_ = leaf_temp;
   atm_kpa_ = atm_kpa;
   atm_o2_kpa_ = atm_o2_kpa;
   PPFD_ = PPFD;
   soil_moist_ = soil_moist;
   root_collar_psi_ = root_collar_psi;
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

// Vulnerability curve for leaf (parameterised from Potkay et al. 2021) [frac]
double VC_l(double psi){
  double b_r = 0.85; //Mpa
  double c_r = 0.81;
  return(exp(-(pow(-psi/b_r,c_r))));
}
// Vulnerability curve for root (parameterised from Potkay et al. 2021) [frac]
double vulnerability_curve_root(double P_soil){
  double b_r = 1.29; //Mpa
  double c_r = 2.65;
  return(exp(-(pow(-P_soil/b_r,c_r))));
}

// Vulnerability curve for sapwood (parameterised from Potkay et al. 2021) [frac]
double VC_sw(double psi){
  double b_r = 5.32; //Mpa
  double c_r = 0.80;
  return(exp(-(pow(-psi/b_r,c_r))));
}

double Leaf::E_from_Soil_to_Root_Collar(double P_x_r, std::vector<double> P_soil, double n_soil, std::vector<double> z_soil_mid, double dz, double LA, std::vector<double> c_r_H, std::vector<double> c_r_V, double beta_R_H, double beta_R_V){
    //   // Integration steps
    // std::cout << "n_soil: " << n_soil << "z_soil_mid: " << z_soil_mid.size() << "P_soil1: " << P_soil[0] << "P_soil2: " << P_soil[1] << "P_soil3: " << P_soil[2]<< "P_soil4: " << P_soil[3] << "P_soil5: " << P_soil[4]<< std::endl;
    // number of intergration steps
    double n = 20;
    // water density
    double water_dens = 1e3;
    // gravity // m s^-2
    double gravity = 9.8;    
    
    // root vulnerability curve parameters
    double b_root = 1.29; //Mpa
    double c_root = 2.65;

    // total root carbon, make it same size as crb
    std::vector<double> c_r;
    c_r.reserve(c_r_V.size());

    // add carbon from horizantal to vertical components
    for (size_t i = 0; i < c_r_V.size(); i++)
    {
      c_r[i] = c_r_H[i] + c_r_V[i];
      // std::cout << "c_r[i]:" << c_r[i] << std::endl;
    }

    // horizantal root resistance
    std::vector<double> r_R_H_min;
    r_R_H_min.reserve(c_r_V.size());
    
    for (size_t i = 0; i < c_r_V.size(); i++)
    {
      r_R_H_min[i] = beta_R_H/c_r_H[i];
    }
    
    // vertical root resistance
    std::vector<double> r_R_V;
    
    // Make sure r_R_V has the correct size
    r_R_V.resize(c_r_V.size());
    
    for (size_t i = 0; i < c_r_V.size(); i++)
    {
      r_R_V[i] = beta_R_V * ((dz * dz) / c_r_V[i]);
      // std::cout << "r_R_V[i]:" << r_R_V[i] << std::endl;
    }
    
    // cumulative vertical sum of root resistance
    std::vector<double> r_R_V_sum;
    
    // Cumulative sum of vertical root resistance
    std::partial_sum(r_R_V.begin(), r_R_V.end(), std::back_inserter(r_R_V_sum));
    std::cout << n_soil << std::endl;
    // Set up vector of root water uptake from layer
    // std::vector<double> soil_consumption_;
    soil_consumption_.reserve(n_soil);
    soil_consumption_.resize(n_soil);

    // Set up vector of fractional resistance from each layer
    std::vector<double> f_r;
    
    for(size_t i = 0; i < n_soil; i++){

      std::cout << i << "inside loop" << std::endl;

    // Find the most negative soil potential out of the given soil layer and the root collar
    double P_src_min = std::min(P_soil[i], P_x_r);
    // Find the least negative soil potential out of the given soil layer and the root collar
    double P_src_max = std::max(P_soil[i], P_x_r);
    
    // std::cout << "P_src_min: " << P_src_min << "P_src_max: "<< P_src_max << std::endl;

     // If root collar soil water potential equals the soil water potential in a given layer
    if(P_x_r == P_soil[i]){
          std::cout << i << "first opinside loop" << std::endl;

      // Fraction of conductance in roots in a given layer at most negative soil water potential
      double f_ri = exp(-(pow(-P_src_min/b_root,c_root)));
      // std::cout << "f_ri, px = ps: "<< f_ri << std::endl;

      // Fraction of conductance in roots in a given layer at most negative soil water potential
      double r_R_H = r_R_H_min[i] / f_ri; // [MPa * s * (mol H2O)^-1]
      // std::cout << "r_R_H:" << r_R_H << "r_R_H_min[i]" << r_R_H_min[i] << "f_ri: " << f_ri << std::endl;

      // Total root resistance (horizantal plus vertical)
      double r_R = r_R_H + r_R_V_sum[i];
      // std::cout << "r_R:" << r_R << "r_R_V_sum[i]" << r_R_V_sum[i] << std::endl;

      // Transpiration is equivalent to gravitational water loss (i.e. layer gains water)
      double E_i = -(water_dens * gravity * z_soil_mid[i] / 1e6) / r_R / LA;
      // std::cout << "E_i:" << E_i << std::endl;

      soil_consumption_.push_back(E_i);
      f_r.push_back(f_ri);

    }
    else if((P_soil[i] - P_x_r) == (water_dens * gravity * z_soil_mid[i] / 1e6)){
                std::cout << i << "second opinside loop" << std::endl;

      // If pressure difference perfectly balances gravity transpiration is equal to zero
      double E_i = 0; // [mol H2O / m^2 / s]
      // std::cout << "E_i:" << E_i << std::endl;

      double f_ri = exp(-(pow(-P_src_min/b_root,c_root)));
      // std::cout << "f_ri, gravity: "<< f_ri << std::endl;

      soil_consumption_[i] = E_i;
      f_r[i] = f_ri;

    } else{
                      std::cout << i << "third opinside loop" << std::endl;

      // Sequence through the most negative to least negative soil water potential
      
      double step = (P_src_max - P_src_min)/n;
      std::vector<double> f_r;
      f_r.reserve(n+1);

      // std::cout << "n:" << n << std::endl;
      for (double P_src_step = P_src_min; P_src_step <= P_src_max; P_src_step += step) {
        // std::cout << "P_src_step: " << P_src_step << std::endl;
       if(P_src_step > 0){
        f_r.push_back(exp(-(pow(0/b_root,c_root))));
      // std::cout << "f_ri, working, 0: "<< exp(-(pow(0/b_root,c_root))) << std::endl;

      } else{
        f_r.push_back(exp(-(pow(-P_src_step/b_root,c_root))));
              // std::cout << "f_r[1]: "<< f_r[i] << std::endl;

      // std::cout << "f_ri, working: "<< exp(-(pow(-P_src_step/b_root,c_root))) << std::endl;
      }
    }
    double f_r_sum = 0;

    // Find the average f_ri

    // std::cout << "size" << f_r.size()<< std::endl;

    for (size_t i = 0; i < n_soil; i++ ){
      // std::cout << "f_r[i]_to_sum" << f_r.at(i) << std::endl;
      f_r_sum += f_r.at(i);
            // std::cout << "f_r_sum" << f_r_sum << std::endl;

    }


    double f_r_average = f_r_sum / n;
    // std::cout << "f_ri_average:" << f_r_average<< "f_r_sum:" << f_r_sum<< "n:" << n << std::endl;

    // Find the horizantal resistance in a given layer by dividing the minimum resistance (i.e. maximum conductivity) by the fractional loss of conductivity
    double r_R_H = r_R_H_min[i] / f_r_average; // [MPa * s * (mol H2O)^-1]
    // std::cout << "r_R_H:" << r_R_H<< "r_R_H_min[i]:" << r_R_H_min[i] << "f_r_average:" << f_r_average << std::endl;

    // Find the total resistance in a given layer by adding the vertical resistance in that layer
    double r_R = r_R_H + r_R_V_sum[i]; // [MPa * s * (mol H2O)^-1]
    // std::cout << "r_R:" << r_R << "r_R_V_sum[i]" << r_R_V_sum[i] << std::endl;

    // Transpiration is equal to the potentail gradient between the root collar and the soil, accounting for gravitational potential
    double E_i = (P_soil[i] - P_x_r - water_dens * gravity * z_soil_mid[i]/ 1e6) / r_R / LA; // [mol H2O / m^2 / s]
    // std::cout << "E_i:" << E_i << std::endl;

    root_collar_psi_ = P_x_r;
    soil_consumption_[i] = E_i;
    f_r[i] = f_r_average;
    }
  }
  
  std::vector<double> r_R_H;
  // Total transpiration equal to sum of uptake from each layer

   E_up_ = 0;

  for(size_t i = 0; i < n_soil; i++){
  std::cout << "soil_consumption_.size():" << soil_consumption_.size() << std::endl;  
  std::cout << "soil_consumption_:" << i << ":" << soil_consumption_[i] << std::endl;
    E_up_ += soil_consumption_[i];
  }
//  std::cout << "E_up_" << E_up_ << std::endl;

  std::vector<double> r_R;
  // Recalculate resistances in each layer (TODO: Bit unsure about why z_soil_mid is [i])
  for (size_t i = 0; i < z_soil_mid.size(); i++){


    r_R.push_back((P_soil[i] - P_x_r - water_dens * gravity * z_soil_mid[i] / 1e6) / soil_consumption_[i] / LA);
    double r_R_H_max = std::max((r_R[i] - r_R_V_sum[i]), r_R_H_min[i]);
    //  std::cout << "P_soil[i]: " << P_soil[i]<< "P_x_r: " << P_x_r << "E_soil[i]: " << E_soil[i] << "z_soil_mid[i]:" << z_soil_mid[i] << "r_R[i]: " << r_R[i] << "r_R_V_sum[i]: " << r_R_V_sum[i] << "r_R_H_min[i]: " << r_R_H_min[i] << std::endl;

    r_R_H.push_back(r_R_H_max);
  }
  return E_up_;
// std::cout << "gravity:" << gravity << std::endl;
}



// A - gc curves

// returns difference between co-limited assimilation and stom_cond_CO2, to be minimised (umol m^-2 s^-1)
double Leaf::E_column(double x, std::vector<double> psi_soil_seperate_, double psi_leaf) {
  std::cout << "starting E column" << std::endl;

  double E_soil_to_root = E_from_Soil_to_Root_Collar(x, psi_soil_seperate_,2, {0.05, 0.15}, 0.1, 1, {20,30}, {20,30}, 3.4e3,  9.4e4);
  // double E_soil_to_root = -0.2;
  root_collar_psi_ = -x;
  // double E_root_to_leaf = x;
  double E_root_to_leaf = transpiration(psi_leaf);
  // std::cout << "root_collar_psi_:" << root_collar_psi_ << "E_soil_to_root: " << E_soil_to_root << "E_root_to_leaf: " <<E_root_to_leaf  << "psi_crit: " << psi_crit << std::endl;

  return E_soil_to_root - E_root_to_leaf;
}

double Leaf::E_column_zero(double x, std::vector<double> psi_soil_seperate_) {

  double E_soil_to_root = E_from_Soil_to_Root_Collar(x, psi_soil_seperate_,2, {0.05, 0.15}, 0.1, 1, {20,30}, {20,30}, 3.4e3,  9.4e4);
    // std::cout << "starting E column:" << E_soil_to_root << std::endl;

  return E_soil_to_root;
}

// converts psi stem to ci, used to find ci which makes A(ci) = gc(ca - ci)
double Leaf::find_root_psi(double wettest_soil_layer, std::vector<double> psi_soil_seperate_, double psi_leaf, int find_root_crit) {
  // not clear what x is here
  

  auto target = [&](double x) mutable -> double {
    if(find_root_crit == 1){
    return E_column(x, psi_soil_seperate_, psi_leaf);
    } else{
      return E_column_zero(x, psi_soil_seperate_);
    }
  };
  // std::cout << "-psi_crit:" << -psi_crit << "wettest_soil_layer:" << wettest_soil_layer << std::endl;

  std::cout << "-psi_crit:" <<-psi_crit << std::endl;

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  return util::uniroot(target, -psi_crit, wettest_soil_layer, 1e-11, ci_niter);

}

double Leaf::find_psi_stem_from_psi_root(double psi_root, std::vector<double> psi_soil_seperate_){
  double E_soil_to_root = E_from_Soil_to_Root_Collar(psi_root, psi_soil_seperate_,2, {0.05, 0.15}, 0.1, 1, {20,30}, {20,30}, 3.4e3,  9.4e4);
  double psi_stem = transpiration_to_psi_stem(E_soil_to_root);
  return psi_stem;
}

double Leaf::find_root_collar_psi(std::vector<double> soil_moist_){
double soil_number_of_depths = soil_moist_.size();
std::vector<double> psi_soil_seperate_;
psi_soil_seperate_.reserve(soil_number_of_depths);
TF24_Environment environment;
// std::cout << environment.soil_number_of_depths << std::endl;
// std::cout << -environment.psi_from_soil_moist(soil_moist_[0]) << std::endl;

  for (int a = 0; a < soil_number_of_depths; a++){
    // std::cout << "a:" << a << std::endl;
    // std::cout << soil_moist_[a] << std::endl;
    psi_soil_seperate_.push_back(-environment.psi_from_soil_moist(soil_moist_[a])/1e6);  
    // std::cout << psi_soil_seperate_[a] << std::endl;
  }
    // double wettest_soil_layer = -0.1;

  double wettest_soil_layer = *std::max_element(psi_soil_seperate_.begin(), psi_soil_seperate_.end());
// std::cout << "wettest_soil_layer:" << wettest_soil_layer << std::endl;
std::cout << "startroot_crit!!!:" << std::endl;

double root_crit = find_root_psi(wettest_soil_layer, psi_soil_seperate_, psi_crit, 1);
std::cout << "root_crit!!!:" << root_crit << std::endl;

double root_zero_E = find_root_psi(wettest_soil_layer, psi_soil_seperate_, psi_crit, 0);
std::cout << "root_zero_E!!!:" << root_zero_E << std::endl;
// std::cout << "transpiration from root zero" << E_from_Soil_to_Root_Collar(root_zero_E, psi_soil_seperate_,2, {0.05, 0.15}, 0.1, 1, {20,30}, {20,30}, 3.4e3,  9.4e4) << std::endl;

double gr = (sqrt(5) + 1) / 2;
  // opt_psi_stem_ = psi_soil_;

  if (-wettest_soil_layer > psi_crit){
    profit_ = 0;
    // profit_ = profit_psi_stem_TF(psi_soil_);
    return profit_;
  }

  // optimise for stem water potential
    double bound_a = -root_zero_E;
    double bound_b = -root_crit;
    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;
    while (abs(bound_b - bound_a) > GSS_tol_abs) {

      std::cout << "bound_c: " << bound_c << std::endl;
      std::cout << "bound_d: " << bound_c << std::endl;
      double psi_stem_c = find_psi_stem_from_psi_root(-bound_c, psi_soil_seperate_);
      root_collar_psi_ = -root_collar_psi_;
      double profit_at_c =
          profit_psi_stem_TF(psi_stem_c);

      double psi_stem_d = find_psi_stem_from_psi_root(-bound_d, psi_soil_seperate_);
      root_collar_psi_ = -root_collar_psi_;
      double profit_at_d =
          profit_psi_stem_TF(psi_stem_d);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    std::cout << "made it out" << std::endl;

    double opt_root_psi = ((bound_b + bound_a) / 2);
    opt_psi_stem_ = find_psi_stem_from_psi_root(-opt_root_psi, psi_soil_seperate_);

    root_collar_psi_ = -root_collar_psi_;
    profit_ = profit_psi_stem_TF(opt_psi_stem_);

    return opt_psi_stem_;
}


// std::vector<double> Leaf::root_collar_psi(std::vector<double> soil_moist_){

// // Density of water
// double rho = 1e3; //kg m^-3

// // Acceleration due to gravity
// double g = 9.8; // m s^-2
// double soil_number_of_depths = soil_moist_.size();
//   std::vector<double> psi_soil_;
// psi_soil_.reserve(soil_number_of_depths);
// TF24_Environment environment;
// std::cout << environment.soil_number_of_depths << std::endl;
// std::cout << -environment.psi_from_soil_moist(soil_moist_[0]) << std::endl;

//   for (int a = 0; a < soil_number_of_depths; a++){
//     std::cout << "a:" << a << std::endl;
//     std::cout << soil_moist_[a] << std::endl;
//     psi_soil_.push_back(-environment.psi_from_soil_moist(soil_moist_[a])/10e6);  
//     std::cout << psi_soil_[a] << std::endl;
//   }

//   std::vector<double>::iterator min_it = std::min_element(psi_soil_.begin(), psi_soil_.end());
//   int min_index = std::distance(psi_soil_.begin(), min_it);
//   double min_index_double = min_index;
//   double depth_min_layer = min_index_double*environment.delta_z;

//   std::cout << "minIndex:" << min_index << std::endl;

//     std::vector<double>::iterator max_it = std::max_element(psi_soil_.begin(), psi_soil_.end());
//   int max_index = std::distance(psi_soil_.begin(), max_it);
//   std::cout << "maxIndex:" << max_index << std::endl;



//   //     Start defining the boundary conditions for the bisection search for the optimum root collar pressure (P_x_r)
//   // The upper value (i.e. least negative value) is equal to the wettest soil layer
  
//   double P_x_r_UB_0 = psi_soil_[max_index];
//   double P_x_r_LB_0 = psi_soil_[min_index] - rho * g * depth_min_layer / 1e6;

//   if(abs(P_x_r_LB_0 - P_x_r_UB_0) < 1){
//     P_x_r_LB_0 = P_x_r_UB_0 - 1;
//   }
//   // Find the middle root collar pressure between the two boundaries
//   double P_x_r_M_0 = (P_x_r_UB_0 + P_x_r_LB_0) / 2;
  
//   // Set up a vector of root collar pressure
//   std::vector<double> P_x_r_3_0;
//   P_x_r_3_0.push_back(P_x_r_UB_0);
//   P_x_r_3_0.push_back(P_x_r_M_0);
//   P_x_r_3_0.push_back(P_x_r_LB_0);
  
//   // Create a vector of transpiration values
//   std::vector<double> E_3_0(3,0);
// std::cout << !(std::any_of(E_3_0.begin(), E_3_0.end(), [](double v){ return v < 0; }) &&
//          std::any_of(E_3_0.begin(), E_3_0.end(), [](double v){ return v > 0; })) << std::endl;
//     //While loop continues until transpiration vector includes a negative and a positive value
//   while (!(std::any_of(E_3_0.begin(), E_3_0.end(), [](double v){ return v < 0; }) &&
//          std::any_of(E_3_0.begin(), E_3_0.end(), [](double v){ return v > 0; }))){
//           for(size_t i = 0; i < P_x_r_3_0.size(); i++){
//                       std::cout << i << std::endl;
//             double P_x_r_ = P_x_r_3_0[i];
//             std::cout << P_x_r_ << std::endl;
//             E_3_0.at(i) = E_from_Soil_to_Root_Collar(P_x_r_, psi_soil_,2, {0.05, 0.15}, 0.1, 1, {20,30}, {20,30}, 3.4e3,  9.4e4);
          
//           std::cout << "E_3_0[0]" << E_3_0[0]<< "E_3_0[1]" << E_3_0[1]<< "E_3_0[2]" << E_3_0[2] << std::endl;

//             //   E_3_0[i] = E_up_;
//           }    // loop body
//           std::cout << "start_adding" << std::endl;

//             // If no negative transpiration values then increase upper boundary root collar psi by an increment
//       int count_lower = std::count_if(E_3_0.begin(), E_3_0.end(), [](double value) {
//         return value < 0;
//     });
//           int count_higher = std::count_if(E_3_0.begin(), E_3_0.end(), [](double value) {
//         return value > 0;
//     });
        
// std::cout << "count_lower" << count_lower << "count_higher" << count_higher << std::endl;
  
// if(count_lower == 0){
//     P_x_r_3_0[0] = P_x_r_3_0[0] + 0.1;
//   }
  
// if(count_higher == 0){
//     P_x_r_3_0[2] = P_x_r_3_0[2] - 0.1;
//   }
//   std::cout << "P_x_r_3_0[0]" << P_x_r_3_0[0] << "P_x_r_3_0[1]" << P_x_r_3_0[1] << "P_x_r_3_0[2]" << P_x_r_3_0[2] << std::endl;
//   // Define new midsection root collar value
//   P_x_r_3_0[1] = (P_x_r_3_0[0] + P_x_r_3_0[2]) / 2;
// }  std::cout << "P_x_r_3_0[0]" << P_x_r_3_0[0] << "P_x_r_3_0[1]" << P_x_r_3_0[1] << "P_x_r_3_0[2]" << P_x_r_3_0[2] << std::endl;

//     return P_x_r_3_0;
// }

  //   Start defining the boundary conditions for the bisection search for the optimum root collar pressure (P_x_r)
  // The upper value (i.e. least negative value) is equal to the wettest soil layer




  // // While loop continues until transpiration vector includes a negative and a positive value
  // while (!(std::any_of(E_3_0.begin(), E_3_0.end(), [](double v){ return v < 0; }) &&
  //        std::any_of(E_3_0.begin(), E_3_0.end(), [](double v){ return v > 0; }))){
  //         for(size_t i = 0; i < P_x_r_3_0.size(); i++){
  //           double P_x_r = P_x_r_3_0[i];
  //           E_from_Soil_to_Root_Collar();
  //           E_3_0[i] <- E_up_;
  //         }    // loop body
  //         }

  // return E_3_0;
  







// while((length(E_3_0[E_3_0 < 0]) == 0) +
//       (length(E_3_0[E_3_0 > 0]) == 0) > 0){
  
//   # iterate across the three root collar values
//   for(i in 1:length(P_x_r_3_0)){
    
//     # start with upper boundary first, through to lower
    
//     P_x_r_i = P_x_r_3_0[i]
    
//     # calculate soil to root transpiration value for a given root collar psi
//     E_3_0[i] <- E_from_Soil_to_Root_Collar(P_x_r_i, P_soil = P_soil, rho = rho, g = g)[[1]]
//   }
  
//   # If no negative transpiration values then increase upper boundary root collar psi by an increment
//   if(length(E_3_0[E_3_0 <0]) == 0){
//     P_x_r_3_0[1] = P_x_r_3_0[1] + 0.1;
//   }
  
//   # If no postive transpiration values then decrease lower boundary root collar psi by an increment
//   if(length(E_3_0[E_3_0 >0]) == 0){
//     P_x_r_3_0[3] = P_x_r_3_0[3] - 0.1;
//   }
  
//   # Define new midsection root collar value
//   P_x_r_3_0[2] = (P_x_r_3_0[1] + P_x_r_3_0[3]) / 2;
// }

//   return psi_soil_;
// }

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
  
  return leaf_specific_conductance_max_ * integrator.integrate(f, root_collar_psi_, psi_stem);
 }

//calculates supply-side transpiration from psi_stem and root_collar_psi_, returns kg h20 s^-1 m^-2 LA
double Leaf::transpiration(double psi_stem) {
  // std::cout << "leaf_specific_conductance_max_: " << leaf_specific_conductance_max_ << "psi_stem: " << psi_stem << "root_collar_psi_: "  << root_collar_psi_<<  "transpiration_from_psi.eval(root_collar_psi_):" << transpiration_from_psi.eval(root_collar_psi_)<< "transpiration_from_psi.eval(psi_stem):"  << transpiration_from_psi.eval(psi_stem) << "root_collar_psi__internal_ " << root_collar_psi_ << std::endl;
  // integration of proportion_of_conductivity over [root_collar_psi_, psi_stem]
  return leaf_specific_conductance_max_ * (transpiration_from_psi.eval(psi_stem) - transpiration_from_psi.eval(root_collar_psi_));
  // return (transpiration_full_integration(psi_stem));

  
}

// converts a known transpiration to its corresponding psi_stem, returns -MPa
double Leaf::transpiration_to_psi_stem(double transpiration_) {
  // integration of proportion_of_conductivity over [root_collar_psi_, psi_stem]

  // std::cout << "root_collar_psi_from_transpiration_to_psi_stem: " << root_collar_psi_ << std::endl;
  // std::cout << "transpiration_from_psi.eval(root_collar_psi_) full:" << transpiration_from_psi.eval(-root_collar_psi_) << std::endl;

  // std::cout << "transpiration_:" << transpiration_ << std::endl;
  // std::cout << "leaf_specific_conductance_max_:" << leaf_specific_conductance_max_ << std::endl;

  double E_psi_stem = transpiration_/leaf_specific_conductance_max_ +  transpiration_from_psi.eval(-root_collar_psi_);

    // std::cout << "E_psi_stem:" << E_psi_stem << std::endl;


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


// A - gc curves

// returns difference between co-limited assimilation and stom_cond_CO2, to be minimised (umol m^-2 s^-1)
double Leaf::assim_minus_stom_cond_CO2(double x, double psi_stem) {

  double assim_colimited_x_ = assim_colimited(x);

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
  
  if (root_collar_psi_ >= psi_stem){
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


// Hydraulic cost equations

// Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::hydraulic_cost_Sperry(double psi_stem) {
  double k_l_soil_ = leaf_specific_conductance_max_ * proportion_of_conductivity(root_collar_psi_);
  double k_l_stem_ = leaf_specific_conductance_max_ * proportion_of_conductivity(psi_stem);
  
  hydraulic_cost_ = k_l_soil_ - k_l_stem_;
  
  return hydraulic_cost_;
}

double Leaf::hydraulic_cost_TF(double psi_stem) {
//hydraulic_cost_ = 1e6 * 
  //  hk_s /(365*24*60*60)* 
    //(1/a_bio_) * 
    //rho_ * sapwood_volume_per_leaf_area_ * pow((1 - proportion_of_conductivity(psi_stem)), beta2);

  hydraulic_cost_ = g1_TF24 * pow((1 - proportion_of_conductivity(psi_stem)), beta2);


return hydraulic_cost_;
}

// Profit functions

double Leaf::profit_psi_stem_Sperry(double psi_stem) {

set_leaf_states_rates_from_psi_stem(psi_stem);

  double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem);

  return benefit_ - lambda_ * hydraulic_cost_;
}


double Leaf::profit_psi_stem_TF(double psi_stem) {
set_leaf_states_rates_from_psi_stem(psi_stem);

  double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_TF(psi_stem);

  return benefit_ - hydraulic_cost_;
}


//optimisation functions


// need docs on Golden Section Search.
void Leaf::optimise_psi_stem_Sperry() {

  double gr = (sqrt(5) + 1) / 2;
  opt_psi_stem_ = root_collar_psi_;


  if ((PPFD_ < 1.5e-8 )| (root_collar_psi_ > psi_crit)){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    return;
  }

  // optimise for stem water potential
    double bound_a = root_collar_psi_;
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
  

void Leaf::optimise_psi_stem_TF() {

  double gr = (sqrt(5) + 1) / 2;
  opt_psi_stem_ = root_collar_psi_;

  if (root_collar_psi_ > psi_crit){
    profit_ = profit_psi_stem_TF(root_collar_psi_);
    return;
  }

  // optimise for stem water potential
    double bound_a = root_collar_psi_;
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
