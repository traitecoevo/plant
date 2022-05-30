// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/canopy.h>
#include <plant/interpolator.h>

using namespace Rcpp;

namespace plant {

// double check best namespace for constants (private vs global)
static const double PAR_to_SW = 1.0; // dummy numbers
static const double slp = 1.0;
static const double gamma = 1.0;
static const double lh = 1.0;
static const double MH2O = 1.0;
static const double g_to_kg = 1000.0;
static const double kg_to_m3 = 1000.0;
static const double sec_2_hr = 3600.0;

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
  double PPFD;
  double T;
  double r_soil;
  double theta_wp;
  double theta_fc;
  double swf;

  virtual void compute_rates(std::vector<double> const& resource_depletion) {
    double infiltration;
    double evaporation;
    double drainage;
    double net_flux;

    // canopy_openness -> shading_above -> includes k_I
    double ground_radiation = PPFD * get_environment_at_height(0);;

    // treat each soil layer as a separate resource pool
    for (size_t i = 0; i < vars.state_size; i++) {

      if(i == 0) {
        infiltration = extrinsic_drivers["rainfall"].eval(time);

        // Evaporation at soil surface
        double E_bare_soil_pot_mol = (1.0 - r_soil) * PPFD * PAR_to_SW * slp / ((slp + gamma) * lh);

        // mols ->  m3/m2/s
        double E_bare_soil_pot_m = std::max(0.0, E_bare_soil_pot_mol * MH2O * g_to_kg * kg_to_m3);

        // need to access total leaf area for environment?
        double patch_total_area_leaf = 1.0;
        double Ev_pot = E_bare_soil_pot_m * std::exp(-0.398 * patch_total_area_leaf);

        double soil_wetness = std::pow(((vars.state(i) - theta_wp) / (theta_fc - theta_wp)), swf);

        evaporation = std::max(0.0, Ev_pot * soil_wetness) * sec_2_hr;
      } else {
        evaporation = 0.0;
      }

      // TODO: add drainage
      // if (i == n_soil_layers) {
      //    drainage = something
      // } else {
      //    drainage = vars.state(i + 1) * something;
      // }

      drainage = 0.0;

      net_flux = infiltration - resource_depletion[i] - evaporation - drainage;

      vars.set_rate(i, net_flux);
    }
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
    const double a_psi = 2539.246;
    const double n_psi = 2.130148;
    const double theta_sat = 0.48;

    double psi = std::pow(a_psi * (theta / theta_sat), -n_psi);

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
