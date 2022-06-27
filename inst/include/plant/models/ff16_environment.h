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
static const double lh = 44002.59;
static const double MH2O = 18.02;
static const double g_to_kg = 1/1000.0;
static const double kg_to_m3 = 1/1000.0;

class FF16_Environment : public Environment {
public:
  // constructor for R interface - default settings can be modified
  // except for soil_number_of_depths and canopy_rescale_usually
  // which are only updated on construction
  FF16_Environment(bool canopy_rescale_usually = false,
                   int soil_number_of_depths = 0)
      : canopy_rescale_usually(canopy_rescale_usually) {
    time = 0.0;
    canopy = Canopy();
    vars = Internals(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, 0.0));
    extrinsic_drivers["rainfall"] = interpolator::Interpolator();
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

  // soil variables
  double PPFD = 1000;

  //paramaterised for sandy loam
  double theta_sat = 0.453;
  double r_soil = 0.01;
  double theta_wp = 0.081;
  double theta_fc = 0.180;
  double swf = 1;
  double k_sat = 440.628;
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
        saturation = std::max(0.0, 1 - std::pow(vars.state(i)/theta_sat, b_infil));
        infiltration = extrinsic_drivers["rainfall"].eval(time) *  saturation;

        // // Evaporation at soil surface
        // double E_bare_soil_pot_mol = (1.0 - r_soil) * PPFD * PAR_to_SW * slp / ((slp + gamma) * lh);

        // // mols ->  m3/m2/s
        // double E_bare_soil_pot_m = std::max(0.0, E_bare_soil_pot_mol * MH2O * g_to_kg * kg_to_m3);

        // // need to access total leaf area for environment?
        // double patch_total_area_leaf = 1.0;
        // double Ev_pot = E_bare_soil_pot_m * std::exp(-0.398 * patch_total_area_leaf);

        // double soil_wetness = std::pow(((vars.state(i) - theta_wp) / (theta_fc - theta_wp)), swf);

        // evaporation = std::max(0.0, Ev_pot * soil_wetness)*60*60*24*365;
      } else {
        infiltration = 0.0;
        // evaporation = 0.0;

      }

      // if (i == soil_number_of_depths){
      //   drainage = k_sat * std::pow(vars.state(i)/theta_sat, 2*calc_n_psi() + 3);
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

      drainage = k_sat * std::pow(vars.state(i)/theta_sat, 2*calc_n_psi() + 3);

      // drainage = 0;


      net_flux = infiltration  - drainage -  resource_depletion[i];

        // std::cout << "time: " << time << "; infiltration: " << infiltration << "; theta: " << vars.state(i) << "resource depletion: " << resource_depletion[i] << std::endl;
 
      
      vars.set_rate(i, net_flux);
    }
  }

  double calc_n_psi() const {
  double n_psi = -((log(1500/33))/(log(theta_wp/theta_fc)));
  return n_psi;
}

  double calc_a_psi() const {
  double n_psi = calc_n_psi();
  double a_psi = 1.5e6 * std::pow(theta_wp, calc_n_psi());
  return a_psi;
}

  double calc_psi(double theta_) const {
  double n_psi = calc_n_psi();
  double a_psi = calc_a_psi();

    double psi = a_psi * std::pow(theta_, -n_psi);
  return psi;
}

  double get_psi_soil() const {
    // soil volumetric water: m3.m-3
    // assume one layer for now - later extend to include layers of variable depth
    double theta = get_soil_water_state()[0];

    // later average over all layers
    // for(i in 1:n_soil_layers)
    //    total = sum(environment.vars.state(i))
    //    theta_soil = total / n_soil_layers

    // hardcode for now; later set in enviornment constructor
    double a_psi = calc_a_psi();
    double n_psi = calc_n_psi();

    double psi = calc_psi(theta);

    // TODO: convert theta to psi
    // calc_apsi(theta_fc, theta_wp)*(theta/theta_sat)^(-calc_npsi(theta_fc, theta_wp))

    return psi;
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

inline Rcpp::NumericMatrix get_state(const FF16_Environment environment) {
  return get_state(environment.canopy);
}


}

#endif
