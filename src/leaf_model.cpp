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
    root_c(2.65), //unitless
    root_b(1.29), //-MPa
    root_psi_crit(1.951651), //-MPa 
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
    g1_TF24(7.5), //cost parameter for TF24 profit model umol m^-2 s^-1
    beta_R_H(3.4e2), //proportionality constant between minimum horizontal (intraleyer) root hydraulic resistance and C_r^-1 in [MPa * s * (mol C) / (mol H2O)]
    beta_R_V(9.4e3) //proportionality constant between minimum vertical (interlayer) root hydraulic resistance and dz^2/C_r in [MPa * (mol C) * s / (mol H2O) / m^2]
   {
      setup_transpiration(100); // arg: num control points for integration
      setup_root_vulnerability(100);
      setup_clean_leaf();
}

Leaf::Leaf(double vcmax_25, double c, double b,
           double psi_crit, // derived from b and c,
           double root_c,
           double root_b,
           double root_psi_crit,
           double beta2, double jmax_25, double hk_s,
           double a, double curv_fact_elec_trans, double curv_fact_colim, 
           double GSS_tol_abs,
           double vulnerability_curve_ncontrol,
           double ci_abs_tol,
           double ci_niter,
           double g1_TF24,
           double beta_R_H,
           double beta_R_V)
    : vcmax_25(vcmax_25), // umol m^-2 s^-1 
    c(c), //unitless
    b(b), //-MPa
    psi_crit(psi_crit), //-MPa 
    root_c(root_c), //unitless
    root_b(root_b), //-MPa
    root_psi_crit(root_psi_crit), //-MPa 
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
    g1_TF24(g1_TF24), //cost parameter for TF24 profit model umol m^-2 s^-1
    beta_R_H(beta_R_H),
    beta_R_V(beta_R_V)
   {
      setup_transpiration(vulnerability_curve_ncontrol); // arg: num control points for integration
      setup_root_vulnerability(vulnerability_curve_ncontrol);
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
  area_leaf_ = NA_REAL;
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
  psi_soil_.clear();
  soil_depth_.clear();
  z_soil_mid_.clear();  // ADD THIS LINE
  use_precomputed_z_soil_mid_ = false;
  c_r_V_.clear(); // carbon per layer dedicated to vertical transport (kg m^-2);
  c_r_H_.clear(); // carbon per layer dedicated to horizantal transport (kg m^-2);
  r_R_H_min.clear(); //minimum horizontal portion of root resistance in each soil-layer in [MPa * s * (mol H2O)^-1 m^-2];
  r_R_V.clear(); // vertical root resitance [MPa * s * (mol H2O)^-1 m^-2];
  r_R_V_sum.clear(); // summed vertical root resistance as depth increase;
  soil_consumption_.clear(); // soil consumption mol  m^-2 s^-1;

  soil_number_of_depths_ = NA_INTEGER;
  max_soil_layer = NA_INTEGER; // number of soil layers with root mass greater than 0;
}

//sets various parameters which are constant for a given node at a given time

void Leaf::set_physiology(double area_leaf, const std::vector<double>& mass_root_prop, double rho, double a_bio, double PPFD, const std::vector<double>& psi_soil, const std::vector<double>& soil_depth, double leaf_specific_conductance_max, double atm_vpd, double ca, double sapwood_volume_per_leaf_area, double leaf_temp, double atm_o2_kpa, double atm_kpa) {
    if (psi_soil.size() != soil_depth.size()) {
    util::stop("soil_depth and psi_soil must have the same number of elements");
  }
  area_leaf_ = area_leaf;
  rho_ = rho;
   a_bio_ = a_bio;
   atm_vpd_ = atm_vpd;
   leaf_temp_ = leaf_temp;
   atm_kpa_ = atm_kpa;
   atm_o2_kpa_ = atm_o2_kpa;
   PPFD_ = PPFD;
   psi_soil_ = psi_soil;
   soil_depth_ = soil_depth;
   soil_number_of_depths_ = soil_depth_.size();
   
   if (!(use_precomputed_z_soil_mid_ &&
         z_soil_mid_.size() == static_cast<size_t>(soil_number_of_depths_))) {
     // Fallback for paths that do not provide environment-precomputed midpoints.
     z_soil_mid_.resize(soil_number_of_depths_);
     for (size_t i = 0; i < soil_number_of_depths_; ++i) {
       if (i == 0) {
         z_soil_mid_[i] = (soil_depth_[i] / 2.0);
       } else {
         z_soil_mid_[i] = ((soil_depth_[i - 1] + soil_depth_[i]) / 2.0);
       }
     }
   }

   use_precomputed_z_soil_mid_ = false;
   
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

   dz_ = soil_depth_.back()/soil_number_of_depths_;


  // find max soil layer as last iteration with mass_root_prop greater than 0
  max_soil_layer = 0;
  for (size_t i = 0; i < soil_number_of_depths_; ++i) {
    if (mass_root_prop[i] != 0) {
      max_soil_layer = i + 1;
    }
  }
  c_r_V_.assign(max_soil_layer, 0.0);
  c_r_H_.assign(max_soil_layer, 0.0);
  r_R_H_min.resize(max_soil_layer);
  r_R_V.resize(max_soil_layer);
  r_R_V_sum.resize(max_soil_layer);

  const double dz_sq = dz_ * dz_;
  double vertical_resistance_sum = 0.0;
  for (size_t i = 0; i < max_soil_layer; ++i) {
    const double root_mass = mass_root_prop[i];
    if (root_mass == 0.0) {
      r_R_H_min[i] = 0.0;
      r_R_V[i] = 0.0;
      r_R_V_sum[i] = vertical_resistance_sum;
      continue;
    }

    const double c_r_v = root_mass / 3.0;
    const double c_r_h = root_mass * 2.0 / 3.0;
    c_r_V_[i] = c_r_v;
    c_r_H_[i] = c_r_h;

    // Set horizantal minimum resistance per soil layer (i.e. reciprocal of maximum conductance).
    r_R_H_min[i] = beta_R_H / c_r_h;
    // The vertical conductivity is likely linearly proportional to the root area projected onto the horizontal plane, hence dz^2.
    r_R_V[i] = beta_R_V * dz_sq / c_r_v;
    vertical_resistance_sum += r_R_V[i];
    r_R_V_sum[i] = vertical_resistance_sum;
  }

  // Set up vector of root water uptake from layer
  soil_consumption_.resize(soil_number_of_depths_, 0.0);

  // Find maximum assimilation assuming ci = ca
  assim_max_ = assim_colimited(ca_);
}

// This function calculates the total transpiration from the soil based on the root collar pressure and the respective soil layer pressures
void Leaf::E_from_Soil_to_Root_Collar(double P_x_r, const std::vector<double>& psi_soil){

    E_up_ = 0;

    for(size_t i = 0; i < max_soil_layer; i++){

    // Find the most negative soil potential out of the given soil layer and the root collar
    double P_src_min = std::min(psi_soil[i], P_x_r);

    // Find the least negative soil potential out of the given soil layer and the root collar
    double P_src_max = std::max(psi_soil[i], P_x_r);
    

    if(P_src_min > P_src_max){
    util::stop("P_src_min must be more negative than P_src_max");
    }


     // If root collar soil water potential equals the soil water potential in a given layer
    if(std::abs(P_x_r - psi_soil[i]) < 1e-20){

      // Fraction of conductance in roots in a given layer at most negative soil water potential (but actually is equal to root collar)
      // root_vuln_from_psi is a pre-built spline of exp(-(|psi|/b_root)^c_root)
      double f_ri = root_vuln_from_psi.eval(-P_src_min);

      // Fraction of conductance in roots in a given layer at most negative soil water potential
      double r_R_H = r_R_H_min[i] / f_ri; // [MPa * s * (mol H2O)^-1]

      // Total root resistance (horizantal plus vertical)
      double r_R = r_R_H + r_R_V_sum[i];

      // Transpiration is equivalent to gravitational water loss (i.e. layer gains water)
      double E_i = -(gravity_head * z_soil_mid_[i]) / r_R / area_leaf_ ;

      soil_consumption_[i] = E_i;
      E_up_ += E_i;

    }
    else if((psi_soil[i] - P_x_r) == (gravity_head * z_soil_mid_[i])){

      // If pressure difference perfectly balances gravity transpiration is equal to zero
      double E_i = 0.0; // [mol H2O / m^2 / s]
      
      soil_consumption_[i] = E_i;

      E_up_ += E_i;

    } else{

      // Sequence through the least negative to most negative soil water potential
      // step will be negative
      double step = (P_src_max - P_src_min)/n;

      double f_r_average = 0;

      for (size_t j = 0; j < (n + 1); j++) {
        double stepper = j;
        double P_src_step = P_src_min + step * stepper;
       if(P_src_step > 0){
        // psi > 0 means above-atmospheric pressure; vulnerability = 1 (no loss)
        f_r_average += 1.0 / (n + 1);
      } else{
        // look up pre-computed root vulnerability spline instead of exp(pow(...))
        f_r_average += root_vuln_from_psi.eval(-P_src_step) / (n + 1);
      }
    }

    // Find the horizantal resistance in a given layer by dividing the minimum resistance (i.e. maximum conductivity) by the fractional loss of conductivity
    double r_R_H = r_R_H_min[i] / f_r_average; // [MPa * s * (mol H2O)^-1]

    // Find the total resistance in a given layer by adding the vertical resistance in that layer
    double r_R = r_R_H + r_R_V_sum[i]; // [MPa * s * (mol H2O)^-1]

    // Transpiration is equal to the potentail gradient between the root collar and the soil, accounting for gravitational potential
    double E_i = (psi_soil[i] - P_x_r - gravity_head * z_soil_mid_[i]) / r_R / area_leaf_; // [mol H2O / m^2 / s]

    soil_consumption_[i] = E_i;
    E_up_ += E_i;

    }
  }
  // convert to kg h20 m-2 s-1 consistent with rest of leaf model and environment TODO: possibly change this
  E_up_ = E_up_*0.018015;
}



// This function is used to find root collar pressure which equilibrates the soil-root-stem water continuuum
double Leaf::E_column(double x, const std::vector<double>& psi_soil, double psi_leaf) {


  E_from_Soil_to_Root_Collar(x, psi_soil);
  root_collar_psi_ = -x;
  double E_root_to_leaf = transpiration(psi_leaf, root_collar_psi_);
  return E_up_ - E_root_to_leaf;
}

// This function is used to find root collar pressure where water form soil is equal to zero
double Leaf::E_column_zero(double x, const std::vector<double>& psi_soil) {

  E_from_Soil_to_Root_Collar(x, psi_soil);

  return E_up_;
}

// find root psi based on required condition, i.e. equilibrated continuum, zero water from soil
double Leaf::find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit) {
  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  if (find_root_crit == 1) {
    auto target = [&](double x) -> double {
      return E_column(x, psi_soil, psi_crit);
    };
    return util::uniroot(target, -psi_crit, wettest_soil_layer, 1e-4, ci_niter);
  }

  auto target = [&](double x) -> double {
    return E_column_zero(x, psi_soil);
  };
  return util::uniroot(target, -psi_crit, wettest_soil_layer, 1e-4, ci_niter);

}

// When root pressure is known, find E from soil, then use E from soil to find psi stem
double Leaf::find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil){
  E_from_Soil_to_Root_Collar(psi_root, psi_soil);

  double psi_stem = transpiration_to_psi_stem(E_up_, psi_root);
  return psi_stem;
}

void Leaf::find_root_collar_psi(){

  
  // Psi soil comes in as positive values but is utilised as negative so need to flip TODO: change thi s around
  psi_soil_inverted_.resize(max_soil_layer);
  double wettest_soil_layer = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < max_soil_layer; ++i) {
    const double psi_inverted = -psi_soil_[i];
    psi_soil_inverted_[i] = psi_inverted;
    wettest_soil_layer = std::max(wettest_soil_layer, psi_inverted);
  }

  // Avoid loop if the wettest psi layer is drier than psi_crit in stem, transpiration not possible and so all variables set to 
  // shut down

  if (-wettest_soil_layer >= psi_crit){

    // profit_ = 0;
    root_collar_psi_ = -psi_crit;
    opt_psi_stem_ = psi_crit;
    profit_ = - R_d_ - hydraulic_cost_TF(psi_crit);
    // return profit_;
    return;
  }

  // Avoid loop if the wettest psi layer is drier than psi_crit in stem, transpiration not possible and so all variables set to 
  // shut down
double root_crit = find_root_psi(wettest_soil_layer, psi_soil_inverted_, 1);

// If root crit would have to be larger than psi crit, also avoid loop as above

    if (-root_crit >= psi_crit){
    // profit_ = 0;
    root_collar_psi_ = root_crit;
    opt_psi_stem_ = psi_crit;
    profit_ = - R_d_ - hydraulic_cost_TF(psi_crit);
       return;

  }

// Find root collar where transpiration from soil is 0
double root_zero_E = find_root_psi(wettest_soil_layer, psi_soil_inverted_, 0);

// If assimilation would be less than 0 even at Ca, also end loop
if(assim_max_ < 0){
    opt_psi_stem_ = root_zero_E;
    root_collar_psi_ = root_zero_E;
    E_from_Soil_to_Root_Collar(root_collar_psi_, psi_soil_inverted_);

    profit_ = - R_d_ - hydraulic_cost_TF(-root_collar_psi_);

        if(std::isnan(profit_)){
          util::stop("Error: profit nan");
    }

    return;
}
// opt_psi_stem_ = psi_soil_;


  // optimise for stem water potential
    double bound_a = -root_zero_E;
    double bound_b = std::max(-root_crit,-root_psi_crit);
    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    double psi_stem_c    = find_psi_stem_from_psi_root(-bound_c, psi_soil_inverted_);
    double profit_at_c   = profit_psi_stem_TF(psi_stem_c, bound_c);

    double psi_stem_d    = find_psi_stem_from_psi_root(-bound_d, psi_soil_inverted_);
    double profit_at_d   = profit_psi_stem_TF(psi_stem_d, bound_d);

while (std::abs(bound_b - bound_a) > GSS_tol_abs) {
  if (profit_at_c > profit_at_d) {
    bound_b    = bound_d;
    bound_d    = bound_c;  
    profit_at_d = profit_at_c;  // reuse
    bound_c    = bound_b - (bound_b - bound_a) / gr;
    psi_stem_c = find_psi_stem_from_psi_root(-bound_c, psi_soil_inverted_);
    profit_at_c = profit_psi_stem_TF(psi_stem_c, bound_c);  // 1 new eval
  } else {
    bound_a    = bound_c;
    bound_c    = bound_d;  
    profit_at_c = profit_at_d;  // reuse
    bound_d    = bound_a + (bound_b - bound_a) / gr;
    psi_stem_d = find_psi_stem_from_psi_root(-bound_d, psi_soil_inverted_);
    profit_at_d = profit_psi_stem_TF(psi_stem_d, bound_d);  // 1 new eval
  }
}


    double opt_root_psi = ((bound_b + bound_a) / 2);

    opt_psi_stem_ = find_psi_stem_from_psi_root(-opt_root_psi, psi_soil_inverted_);


    root_collar_psi_ = opt_root_psi;
    profit_ = profit_psi_stem_TF(opt_psi_stem_, root_collar_psi_);

    if(std::isnan(profit_)){
          util::stop("Error: nan");
    }
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

// pre-compute root vulnerability curve f(psi) = exp(-(|psi|/b_root)^c_root) as a spline,
// evaluated over the range [0, psi_max_root] where conductivity drops to 1%.
// This avoids repeated exp(pow(...)) calls inside E_from_Soil_to_Root_Collar.
void Leaf::setup_root_vulnerability(double resolution) {
  auto x_psi_root = std::vector<double>{0.0};
  auto y_f_r       = std::vector<double>{1.0}; // f(0) = exp(0) = 1
  // upper limit: psi where conductivity = 1%
  double psi_max_root = root_b * pow(log(1.0 / 0.01), 1.0 / root_c);
  double step = psi_max_root / resolution;
  for (double psi = step; psi <= psi_max_root; psi += step) {
    x_psi_root.push_back(psi);
    y_f_r.push_back(exp(-pow(psi/root_b, root_c)));
  }
  root_vuln_from_psi.init(x_psi_root, y_f_r);
  root_vuln_from_psi.set_extrapolate(true); // clamp to last value beyond range
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

double Leaf::transpiration_full_integration(double psi_stem, double psi_upstream) {
  std::function<double(double)> f;
  f = [&](double psi) -> double { return proportion_of_conductivity(psi); };
  
  return leaf_specific_conductance_max_ * integrator.integrate(f, psi_upstream, psi_stem);
 }

//calculates supply-side transpiration from psi_stem and root_collar_psi_, returns kg h20 s^-1 m^-2 LA
double Leaf::transpiration(double psi_stem, double psi_upstream) {

  
  // integration of proportion_of_conductivity over [root_collar_psi_, psi_stem]
  return leaf_specific_conductance_max_ * (transpiration_from_psi.eval(psi_stem) - transpiration_from_psi.eval(psi_upstream));
  // return (transpiration_full_integration(psi_stem));

  
}

// converts a known transpiration to its corresponding psi_stem, returns -MPa
double Leaf::transpiration_to_psi_stem(double transpiration_, double psi_upstream) {
  // integration of proportion_of_conductivity over [root_collar_psi_, psi_stem]


  double E_psi_stem = transpiration_/leaf_specific_conductance_max_ +  transpiration_from_psi.eval(-psi_upstream);


  return psi_from_transpiration.eval(E_psi_stem);
  }

// returns stomatal conductance to CO2, mol C m^-2 LA s^-1
double Leaf:: stom_cond_CO2(double psi_stem, double psi_upstream) {
  double transpiration_ = transpiration(psi_stem, psi_upstream);
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
double Leaf::assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream) {

  double assim_colimited_x_ = assim_colimited(x);

  double stom_cond_CO2_x_ = stom_cond_CO2(psi_stem, psi_upstream);
  return assim_colimited_x_ * umol_to_mol -
         (stom_cond_CO2_x_ * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
}

// converts psi stem to ci, used to find ci which makes A(ci) = gc(ca - ci)
double Leaf::psi_stem_to_ci(double psi_stem, double psi_upstream) {
  const double stom_cond_CO2_fixed = stom_cond_CO2(psi_stem, psi_upstream);

  auto target = [&](double x) mutable -> double {
    const double assim_colimited_x_ = assim_colimited(x);
    return assim_colimited_x_ * umol_to_mol -
      (stom_cond_CO2_fixed * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
  };

  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  return ci_ = util::uniroot(target, gamma_ * umol_per_mol_to_Pa, ca_, 1e-7, ci_niter);
}

// given psi_stem, find assimilation, transpiration and stomal conductance to c02
void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream) {

  if (psi_upstream >= psi_stem){
    ci_ = gamma_*umol_per_mol_to_Pa;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    } else{
      if(assim_max_ < 0){
        ci_ = gamma_*umol_per_mol_to_Pa;
        transpiration_ = 0;
        stom_cond_CO2_ = 0;
        } else{
      ci_ = psi_stem_to_ci(psi_stem, psi_upstream);
      transpiration_ = transpiration(psi_stem, psi_upstream);
      stom_cond_CO2_ = atm_kpa_ * transpiration_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
      }
    }
  assim_colimited_ = assim_colimited(ci_);
}


// Hydraulic cost equations

// Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::hydraulic_cost_Sperry(double psi_stem, double psi_upstream) {
  double k_l_soil_ = leaf_specific_conductance_max_ * proportion_of_conductivity(psi_upstream);
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

double Leaf::profit_psi_stem_Sperry(double psi_stem, double psi_upstream) {

set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_Sperry(psi_stem, psi_upstream);

  return benefit_ - lambda_ * hydraulic_cost_;
}


double Leaf::profit_psi_stem_TF(double psi_stem, double psi_upstream) {
set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

double benefit_ = assim_colimited_;
  double hydraulic_cost_ = hydraulic_cost_TF(psi_stem);

  return benefit_ - hydraulic_cost_;
}


//optimisation functions


// need docs on Golden Section Search.
void Leaf::optimise_psi_stem_Sperry() {

    if (!(psi_soil_.size() == 1)) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }

  opt_psi_stem_ = psi_soil_[0];


  if ((PPFD_ < 1.5e-8 )| (psi_soil_[0] > psi_crit)){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    return;
  }

  // optimise for stem water potential
    double bound_a = psi_soil_[0];
    double bound_b = psi_crit;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > GSS_tol_abs) {

      double profit_at_c =
          profit_psi_stem_Sperry(bound_c, psi_soil_[0]);

      double profit_at_d =
          profit_psi_stem_Sperry(bound_d, psi_soil_[0]);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem_ = ((bound_b + bound_a) / 2);
    profit_ = profit_psi_stem_Sperry(opt_psi_stem_, psi_soil_[0]);

  }
  

void Leaf::optimise_psi_stem_TF() {

  if (!(psi_soil_.size() == 1)) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }

  opt_psi_stem_ = psi_soil_[0];

  if (psi_soil_[0] > psi_crit){
    profit_ = profit_psi_stem_TF(psi_soil_[0], psi_soil_[0]);
    return;
  }

  // optimise for stem water potential
    double bound_a = psi_soil_[0];
    double bound_b = psi_crit;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;
    while (abs(bound_b - bound_a) > GSS_tol_abs) {

      double profit_at_c =
          profit_psi_stem_TF(bound_c, psi_soil_[0]);

      double profit_at_d =
          profit_psi_stem_TF(bound_d, psi_soil_[0]);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem_ = ((bound_b + bound_a) / 2);
    profit_ = profit_psi_stem_TF(opt_psi_stem_, psi_soil_[0]);

    return;
  }

} // namespace plant
