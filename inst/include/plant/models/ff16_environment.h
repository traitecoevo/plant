// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/canopy.h>
#include <plant/interpolator.h>

using namespace Rcpp;

namespace plant {

// double check best namespace for constants (private vs global)
static const double PAR_to_SW = 0.4376; // slp, gamma and lh are temp-dependent, but assumed at 25 deg C for now.
static const double slp = 0.1887;
static const double gamma = 0.0674;
static const double lh = 44002.59; // latent heat of vaporisation of water
static const double MH2O = 18.02; // molar mass of water (g mol^-1 H20)
static const double g_to_kg = 1/1000.0; // convert grams to kilograms
static const double kg_to_m3 = 1/1000.0; // converget kilograms of water to volume

class FF16_Environment : public Environment {
public:
  // constructor for R interface - default settings can be modified
  // except for soil_number_of_depths and canopy_rescale_usually
  // which are only updated on construction
  FF16_Environment(bool canopy_rescale_usually = false,
                   int soil_number_of_depths = 0,
                   double PPFD_ = 0)
      : canopy_rescale_usually(canopy_rescale_usually) {
    time = 0.0;
    canopy = Canopy();
    vars = Internals(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, 0.0));
  };


  // Light interface
  bool canopy_rescale_usually;

  // private?
  Canopy canopy;

  // Should this be here or in canopy?
  void set_fixed_environment(double value, double height_max) {
    canopy.set_fixed_canopy(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return canopy.get_canopy_at_height(height);
  }

  double canopy_openness(double height) const {
    return canopy.canopy_openness(height);
  }

  void r_init_interpolators(const std::vector<double>& state) {
    canopy.r_init_interpolators(state);
  }

  
  double PPFD = 1000;
  
  // soil variables - paramaterised for sandy loam
  double soil_moist_sat = 0.453; // saturated soil moisture content (m3 water m^-3 soil) 
  double r_soil = 0.01;
  double soil_moist_wp = 0.081; // soil moisture at wilting point (m3 water m^-3 soil) 
  double soil_moist_fc = 0.180; // soil moisture at field capacity (m3 water m^-3 soil) 
  double swf = 1;
  double k_sat = 440.628; //saturated hydraulic conductivity of soil
  double b_infil = 8;

  virtual void compute_rates(std::vector<double> const& resource_depletion) {
    double saturation;
    double infiltration;
    double evaporation;
    double drainage;
    double net_flux;

    //debugging

    // canopy_openness -> shading_above -> includes k_I
    // double ground_radiation = PPFD * get_environment_at_height(0);

    // treat each soil layer as a separate resource pool
    for (size_t i = 0; i < vars.state_size; i++) {

      if(i == 0) {

        // double k_I = 0.5;

        //     double canopy_openness_at_ground = get_environment_at_height(0);
        //     double total_LAI = log(canopy_openness_at_ground) / (-k_I);
        //     // std::cout << "total_LAI" << total_LAI;
        //     double ground_radiation = PPFD * get_environment_at_height(0);



        saturation = std::max(0.0, 1 - std::pow(vars.state(i)/soil_moist_sat, b_infil));
        infiltration = extrinsic_drivers.evaluate("rainfall", time) * saturation;


        // drainage = 200*18.02/1000000 * pow((0.35*10e-3/psi_from_soil_moist(vars.state(i))), 5)

        // // Evaporation at soil surface
        // double E_bare_soil_pot_mol = (1.0 - r_soil) * ground_radiation * PAR_to_SW * slp / ((slp + gamma) * lh);

        // // mols ->  m3/m2/s
        // double E_bare_soil_pot_m = std::max(0.0, E_bare_soil_pot_mol * MH2O * g_to_kg * kg_to_m3);

        // // need to access total leaf area for environment?
        // double patch_total_area_leaf = 1.0;
        // double Ev_pot = E_bare_soil_pot_m * std::exp(-0.398 * total_LAI);

        // double soil_wetness = std::pow(((vars.state(i) - soil_moist_wp) / (soil_moist_fc - soil_moist_wp)), swf);

        // evaporation = std::max(0.0, Ev_pot * soil_wetness)*60*60*24*365;
      } else {
        infiltration = 0.0;
        
        
        // evaporation = 0.0;

      }

      // if (i == soil_number_of_depths){
      //   drainage = k_sat * std::pow(vars.state(i)/soil_moist_sat, 2*n_psi() + 3);
      // }

      // TODO: add drainage
      // if (i == n_soil_layers) {
      //    drainage = something
      // } else {
      //    drainage = vars.state(i + 1) * something;
      // }

      // drainage = 0.4;

      // net_flux = infiltration - resource_depletion[i] - evaporation - drainage;
      // net_flux = infiltration  - resource_depletion[i] - drainage;
      // net_flux = infiltration   - drainage;

      // drainage = k_sat * std::pow(vars.state(i)/soil_moist_sat, 2*n_psi() + 3);

      // drainage = 0;

      net_flux = infiltration  -  resource_depletion[i];
      // std::cout << "resource_depletion" << resource_depletion[i] << std::endl;

      // net_flux = infiltration - evaporation - drainage -  resource_depletion[i];

        // std::cout << "; soil_moist: " << vars.state(i) << "resource depletion: " << resource_depletion[i] << std::endl;

      // std::cout << "net_flux" << net_flux << std::endl;


      vars.set_rate(i, net_flux);
    }
  }
// calculate n_psi (parameter to convert soil mositure to soil water potential) from soil field capacity and soil wilting point
  double n_psi() const {
  double n_psi_ = -((log(1500/33))/(log(soil_moist_wp/soil_moist_fc)));
  return n_psi_;
}

// calculate a_psi (parameter to convert soil mositure to soil water potential) from n_psi and soil wilting point
  double a_psi() const {
  double n_psi_ = n_psi();
  double a_psi_ = 1.5e6 * std::pow(soil_moist_wp, n_psi_);
  return a_psi_;
}

// convert soil moisture to soil water potential
  double psi_from_soil_moist(double soil_moist_) const {
  double n_psi_ = n_psi();
  double a_psi_ = a_psi();

    double psi = a_psi_ * std::pow(soil_moist_, -n_psi_);
  return psi;
}

// convert soil water potential to soil moisture
  double soil_moist_from_psi(double psi_soil_) const {
  double n_psi_ = n_psi();
  double a_psi_ = a_psi();  
  
  return pow((psi_soil_/a_psi_), (-1/n_psi_));
}


  double get_psi_soil() const {
    // soil volumetric water: m3.m-3
    // assume one layer for now - later extend to include layers of variable depth
    double soil_moist_ = get_soil_water_state()[0];
// std::cout << "soil_moist" << soil_moist << std::endl;
    // later average over all layers
    // for(i in 1:n_soil_layers)
    //    total = sum(environment.vars.state(i))
    //    soil_moist_soil = total / n_soil_layers

    // hardcode for now; later set in enviornment constructor
    double a_psi_ = a_psi();
    double n_psi_ = n_psi();

    double psi_ = psi_from_soil_moist(soil_moist_);
    // TODO: convert soil_moist to psi
    // calc_apsi(soil_moist_fc, soil_moist_wp)*(soil_moist/soil_moist_sat)^(-calc_npsi(soil_moist_fc, soil_moist_wp))

    return psi_;
  }

double get_vpd() const {
    return extrinsic_drivers.evaluate("vpd", time);
  }

double get_co2() const {
    return extrinsic_drivers.evaluate("co2", time);
  }


  std::vector<double> get_soil_water_state() const {
    return vars.states;
  }

  // I wonder if this needs a better name? See also environment.h
  Internals r_internals() const { return vars; }

  // R interface
  void set_soil_water_state(std::vector<double> state) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.set_state(i, state[i]);
    }
  }

  // Core functions
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max) {
    canopy.compute_canopy(f_compute_competition, height_max);
  }

  template <typename Function>
  void rescale_environment(Function f_compute_competition, double height_max) {
    canopy.rescale_canopy(f_compute_competition, height_max);
  }

  void clear_environment() {
    canopy.clear();
  }
};

//inline Rcpp::NumericMatrix get_state(const FF16_Environment environment) {
//  return get_state(environment.canopy);
//}
inline Rcpp::List get_state(const FF16_Environment environment, double time) {
  auto ret = get_state(environment.extrinsic_drivers, time);
  ret["canopy"] = get_state(environment.canopy); // does a full copy of ret, not efficient
  auto const& soil_moist_list = environment.get_soil_water_state();
  auto rcpp_soil_moist_vec = Rcpp::NumericVector(soil_moist_list.begin(), soil_moist_list.end());
  ret["soil_moist"] = rcpp_soil_moist_vec;
  return ret;
}

}

#endif
