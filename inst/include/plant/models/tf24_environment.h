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
  // except for soil_number_of_depths
  // which are only updated on construction
  
  TF24_Environment(bool light_availability_spline_rescale_usually = true,
                   int soil_number_of_depths = 5, 
                   double delta_z = 9999, // not using this
                   double soil_moist_sat = 0.453, // saturated soil moisture content (m3 water m^-3 soil) 
                   double K_sat = 440.628, //saturated hydraulic conductivity of soil
                   double a_psi = 8.7, // not currently being used
                   double n_psi = 4.8, // not currently being used
                   double a_infil = 1, // infiltration switch (0-1), 0 no runoff, 1 runoff
                   double b_infil = 8, // unitless, determines infiltration rate
                   double depth = 1.5)  // total depth of soil (m)
      : delta_z(delta_z),
      soil_moist_sat(soil_moist_sat),
      a_psi(a_psi),
      n_psi(n_psi),
      K_sat(K_sat),
      a_infil(a_infil),
      b_infil(b_infil),
      depth(depth)
  {
    time = 0.0;



    // Shading defaults have lower tolerance which are overwritten for speed
    light_availability = ResourceSpline(
                   1e-4,  // light_availability_spline_tol,
                   17,    // light_availability_spline_nbase,
                   16,    // light_availability_spline_max_depth,
                   true //light_availability_spline_rescale_usually)
                  );

    ExtrinsicDrivers extrinsic_drivers;

    extrinsic_drivers_set_constant("PPFD",1800);
    extrinsic_drivers_set_constant("rainfall",1);
    extrinsic_drivers_set_constant("atm_vpd",1);
    extrinsic_drivers_set_constant("ca",40);
    extrinsic_drivers_set_constant("leaf_temp",25);
    extrinsic_drivers_set_constant("atm_o2_kpa",21);
    extrinsic_drivers_set_constant("atm_kpa",100.5);

    set_soil_number_of_depths(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, soil_moist_sat*0.5));
  };
  
  // Number of cumulative auxilliary variables to track in soil moisture model
  double aux_num = 3;
  
  // Setup soil water distribtuion
  void set_soil_number_of_depths(int n) {
    soil_number_of_depths = n;
    
    vars = Internals(soil_number_of_depths + aux_num);

    z.resize(soil_number_of_depths);
    dz.resize(soil_number_of_depths);
    // positive downwards
    water_flux.resize(soil_number_of_depths);

    delta_z = depth / soil_number_of_depths;

    for (int i = 0; i < soil_number_of_depths; i++)
    {
      z[i] = (i + 1) * delta_z;
    }

    for (int i = 0; i < soil_number_of_depths; i++)
    {
      dz[i] = delta_z;
    }
  }
  int get_soil_number_of_depths() const {return soil_number_of_depths;}

  // TODO: should we use auxilliary in internals
  std::vector<double> water_flux;
  std::vector<double> z;
  std::vector<double> dz;

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline light_availability;

  // Light interface
  bool canopy_rescale_usually;
  //distance between layers
  int soil_number_of_depths;
  double delta_z;

  double depth;
  //saturated soil moisture
  double soil_moist_sat;
  //Saturated soil hydraulic conductivity
  double K_sat;
  double a_psi;
  double n_psi;
  double a_infil;
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
  
  virtual void compute_rates(std::vector<double> const &resource_depletion)
  {
    double water_input;
    double rainfall = extrinsic_drivers.evaluate("rainfall", time);
    double infiltration = rainfall*std::max(0.0, 1 - a_infil*std::pow(vars.state(0)/soil_moist_sat, b_infil));

    // treat each soil layer as a separate resource pool
    for (size_t i = 0; i < soil_number_of_depths; i++)
    {

      // initial representation of drainage; to be improved
      if (i == 0)
      {
        water_input = infiltration;
      }
      else
      {
        // m3 m^-2
        water_input = water_flux[i-1];
      }
        // TODO: m3 m^-2
      water_flux[i] = K_sat*pow(vars.state(i)/soil_moist_sat, 2);
      // this function does runoff

      vars.set_rate(i, (water_input - water_flux[i])/dz[i]); 
    }
      vars.set_rate(soil_number_of_depths, rainfall);
      vars.set_rate(soil_number_of_depths + 1, infiltration);
      vars.set_rate(soil_number_of_depths + 2, water_flux[soil_number_of_depths - 1]);
  }

  // calculate K from K_sat based on theta
  double soil_K_from_soil_theta(double theta) {
    //Eq. 5 Zeng and Decker (2009), ref Clapp and Hornberger (1978)
  return K_sat * std::pow(theta/soil_moist_sat, 2*n_psi + 3);
  }


  // convert soil moisture to soil water potential
  double psi_from_soil_moist(double soil_moist_) const {
    return a_psi * std::pow(soil_moist_/soil_moist_sat, -n_psi);
  }

  // convert soil water potential to soil moisture
  double soil_moist_from_psi(double psi_soil_) const {
    return pow((psi_soil_/a_psi), (-1/n_psi))*soil_moist_sat;
  }

  // Easy wrappers. Cn also use `extrinsic_drivers_evaluate("PPFD", time)

  double get_PPFD()      const { return extrinsic_drivers.evaluate("PPFD", time); }
  double get_atm_vpd()   const { return extrinsic_drivers.evaluate("atm_vpd", time); }
  double get_ca()        const { return extrinsic_drivers.evaluate("ca", time); }
  double get_leaf_temp() const { return extrinsic_drivers.evaluate("leaf_temp", time); }
  double get_atm_o2_kpa()const { return extrinsic_drivers.evaluate("atm_o2_kpa", time); } 
  double get_atm_kpa()   const { return extrinsic_drivers.evaluate("atm_kpa", time); } 


  std::vector<double> get_soil_water_state() const { return {vars.states.begin(), vars.states.end() - aux_num}; }
  std::vector<double> get_soil_water_state_cumulative_flux() const { return {vars.states.end()-aux_num, vars.states.end()}; }
  std::vector<double> get_soil_depths() const { return z; }

  // TODO: I wonder if this needs a better name? See also environment.h
  Internals r_internals() const { return vars; }

  // R interface
  void set_soil_water_state(std::vector<double> state) {
    if(state.size() != (vars.state_size- aux_num)) {
      throw std::invalid_argument("Input vector size does not match soil state size.");
    }
    for (size_t i = 0; i < (vars.state_size); i++) {
      if(i < soil_number_of_depths){
        vars.set_state(i, state[i]);
      } else {
        vars.set_state(i, 0);
      }
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

  virtual Rcpp::List r_get_state() const
  {
    
    // Surely an easier way?
    auto const &soil_depth_list = get_soil_depths();
    auto rcpp_soil_depth_vec = Rcpp::NumericVector(soil_depth_list.begin(), soil_depth_list.end());

    auto const &soil_moist_list = get_soil_water_state();
    auto rcpp_soil_moist_vec = Rcpp::NumericVector(soil_moist_list.begin(), soil_moist_list.end());

    auto const &soil_moist_cumulative_flux_list = get_soil_water_state_cumulative_flux();
    auto rcpp_soil_moist_vec_cumulative_flux = Rcpp::NumericVector(soil_moist_cumulative_flux_list.begin(), soil_moist_cumulative_flux_list.end());

    return Rcpp::List::create(
        // auto ret = get_state(environment.extrinsic_drivers, time);

        _["light_availability"] = light_availability.r_get_state(),
        _["soil_moist"] = rcpp_soil_moist_vec,
        _["soil_depth"] = rcpp_soil_depth_vec,
        _["soil_moist_cumulative_flux"] = rcpp_soil_moist_vec_cumulative_flux
    );
  }
  };
}

#endif
