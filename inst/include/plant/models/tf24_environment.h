// Built from  inst/include/plant/models/ff16_environment.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_ENVIRONMENT_H_
#define PLANT_PLANT_TF24_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <plant/interpolator.h>

using namespace Rcpp;

namespace plant {

class TF24_Environment : public Environment {
public:
  // constructor for R interface - default settings can be modified
  // except for soil_number_of_depths and light_availability_spline_rescale_usually
  // which are only updated on construction
  TF24_Environment(bool light_availability_spline_rescale_usually = false,
                   int soil_number_of_depths = 0, 
                   double delta_z = 0.1,
                   double soil_moist_sat = 0.453, // saturated soil moisture content (m3 water m^-3 soil) 
                   double K_sat = 440.628, //saturated hydraulic conductivity of soil
                   double a_psi = 8.7,
                   double n_psi = 4.8,
                   double b_infil = 8)
      : canopy_rescale_usually(canopy_rescale_usually),
      delta_z(delta_z),
      soil_moist_sat(soil_moist_sat),
      a_psi(a_psi),
      n_psi(n_psi),
      K_sat(K_sat),
      b_infil(b_infil) {
    time = 0.0;
    
    light_availability = ResourceSpline();
    light_availability.spline_rescale_usually = light_availability_spline_rescale_usually;

    vars = Internals(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, 0.0));
  };

  //TODO: should we use auxilliary in internals
  std::vector<double> q;

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline light_availability;

  // Light interface
  bool canopy_rescale_usually;
  //distance between layers
  double delta_z;
  //saturated soil moisture
  double soil_moist_sat;
  //Saturated soil hydraulic conductivity
  double K_sat;
  double a_psi;
  double n_psi;
  double b_infil;

  // Ability to prescribe a fixed value
  // TODO: add setting to set other variables like water
  void set_fixed_environment(double value, double height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return light_availability.get_value_at_height(height);
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }

  double PPFD = 1800;

  virtual void compute_rates(std::vector<double> const& resource_depletion) {

    //TODO: add variable depths
    //TODO: track inflow and outflow for mass convservation
    //TODO: check equation for resource consumption

    //Zeng and Decker (2009; https://doi.org/10.1175/2008JHM1011.1), Ireson et al. 2023 (https://doi.org/10.5194/gmd-16-659-2023), can also see Wang et al. 2011 for further agreement
    //Zeng and Decker (2009) suggest the idea of a modified Richards equation (implemented in CABLE in Decker (2015);doi/10.1002/2015MS000507)
    //to account for inadqueate boundary condition in water table but we opt here for simpler version as in Wang et al. (2011; doi/full/10.1029/2010JG001385)
    //as per Ireson et al. (2023) we track soil moisture at a series of nodes and fluxes at the midpoints evenly spaced between the nodes (Figure 1)
    // unlike Ireson et al. (2023) we track theta instead of psi (soil water potential) for conceptual simplicity and to facilitate averaging of layers for transpiration

    
    //number of nodes in soil water column
    int n = vars.state_size;
    //distance between nodes
    // double delta_z = 0.1;
    //helpers
    double dq_dz;
    double dtheta_dt;
    double dpsi_dz;

    //flux across top boundary layer
    q[0] = // rainfall at time
        extrinsic_drivers.evaluate("rainfall", time) * 
        // fraction of precipitation infiltrating
        std::max(0.0, 1 - std::pow(vars.state(0)/soil_moist_sat, b_infil));;
    
    double dtheta_dt_sum = 0;
    double outflow = 0; 
    double resource_depletion_sum = 0;

    //for each node below top node including last node 
    for (size_t i = 0; i < n; i++) {
      //calculate flux across boundary below node (hence i + 1, rather than i)
      if(i < (n-1)){
        //Eq. 3 + 11, Ireson et al. (2023)
        dpsi_dz = (psi_from_soil_moist(vars.state(i + 1)) - psi_from_soil_moist(vars.state(i)))/delta_z;
        q[i + 1] =  0.5 * (soil_K_from_soil_theta(vars.state(i)) + soil_K_from_soil_theta(vars.state(i+1))) * (dpsi_dz - 1);

      } else {
        //Eq. 13 free drainage boundary used at lower boundary
        q[i + 1] = soil_K_from_soil_theta(vars.state(n-1));
        outflow = q[i + 1];
      }
      // Eq. 10, derivative of flux at node w.r.t to depth
      dq_dz = (q[i+1] - q[i])/delta_z;

      // Eq. 2, but we have subtracted resource depletion rates
      dtheta_dt = -dq_dz -  resource_depletion[i];

      dtheta_dt_sum += dtheta_dt;
      resource_depletion_sum += resource_depletion[i];

      vars.set_rate(i, dtheta_dt);
      }

      std::cout << dtheta_dt_sum << "q0" << q[0] << "outflow" << outflow << "resource" << resource_depletion_sum <<  std::endl;

  }

// calculate K from K_sat based on theta
double soil_K_from_soil_theta(double theta) {
  //Eq. 5 Zeng and Decker (2009), ref Clapp and Hornberger (1978)
double K = K_sat * std::pow(theta/soil_moist_sat, 2*n_psi + 3);
return K;
}


// convert soil moisture to soil water potential
  double psi_from_soil_moist(double soil_moist_) const {
    double psi = a_psi * std::pow(soil_moist_/soil_moist_sat, -n_psi);
  return psi;
}

// convert soil water potential to soil moisture
  double soil_moist_from_psi(double psi_soil_) const {
  return pow((psi_soil_/a_psi), (-1/n_psi))*soil_moist_sat;
}

double get_atm_vpd() const {
    return extrinsic_drivers.evaluate("atm_vpd", time);
  }

double get_ca() const {
    return extrinsic_drivers.evaluate("ca", time);
  }

double get_leaf_temp() const {
    return extrinsic_drivers.evaluate("leaf_temp", time);
  }

double get_atm_o2_kpa() const {
    return extrinsic_drivers.evaluate("atm_o2_kpa", time);
  } 

double get_atm_kpa() const {
    return extrinsic_drivers.evaluate("atm_kpa", time);
  } 


  std::vector<double> get_soil_water_state() const {
    return vars.states;
  }

  // TODO: I wonder if this needs a better name? See also environment.h
  Internals r_internals() const { return vars; }

  // R interface
  void set_soil_water_state(std::vector<double> state) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.set_state(i, state[i]);
    }
  }

  // Pre-compute resources available in the environment, as a function of height
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of light_availability spline
    // Note: extinction coefficient was already applied in strategy, so
    // f_compute_competition gives sum of projected leaf area (k L) across species. Just need to apply Beer's law, E = exp(- (k L))
    auto f_light_availability = [&](double height) -> double
    { return exp(-f_compute_competition(height)); };

    // Calculates the light_availability spline, by fitting to the function
    // `f_compute_competition` as a function of height
    light_availability.compute_environment(f_light_availability, height_max, rescale);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }
};


inline Rcpp::List get_state(const TF24_Environment environment, double time) {
  auto ret = get_state(environment.extrinsic_drivers, time);
  ret["light_availability"] = get_state(environment.light_availability);
  return ret;
}
}

#endif
