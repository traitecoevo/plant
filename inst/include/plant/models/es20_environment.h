// Built from  inst/include/plant/models/ff16_environment.h on Wed Sep 16 07:11:11 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_ES20_ENVIRONMENT_H_
#define PLANT_PLANT_ES20_ENVIRONMENT_H_

#include <plant/control.h>
#include <plant/disturbance.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/environment.h>
#include <plant/util.h>
#include <Rcpp.h>

#include <random>
#include <iostream>

using namespace Rcpp;

namespace plant {

class ES20_Environment : public Environment {
public:

  ES20_Environment() {
    // Define an anonymous function to pass got the environment generator
    time = NA_REAL;
    disturbance_regime = 0;
    seed_rain = { 1.0, 1.0, 1.0 };
    seed_rain_index = 0;
    k_I = NA_REAL;
    environment_generator = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
    environment_interpolator = environment_generator.construct(
      [&](double height) {
          return get_environment_at_height(height);
      }, 0, 1); // these are update with init(x, y) when patch is created
    
    k_s = 50;
    stress_mean = 0.75;
    stress_sd = 0.0684931506849315;
    // 
    prepare_stress();

    
  };

  ES20_Environment(double disturbance_mean_interval,
                   std::vector<double> seed_rain_,
                   double k_I_,
                   Control control) {
    k_I = k_I_;
    time = 0.0;
    disturbance_regime = disturbance_mean_interval;
    seed_rain = seed_rain_;
    seed_rain_index = 0;
    environment_generator = interpolator::AdaptiveInterpolator(
      control.environment_light_tol,
      control.environment_light_tol,
      control.environment_light_nbase,
      control.environment_light_max_depth
    );
    environment_interpolator = environment_generator.construct(
      [&](double height) {
          return get_environment_at_height(height);
      }, 0, 1); // these are update with init(x, y) when patch is created
    
    
    //PREPARE STRESS
    k_s = 50;
    stress_mean = control.stress_mean;
    stress_sd = control.stress_sd;
    // 
    k_s = 50;
    stress_mean = 0.75;
    stress_sd = 0.0684931506849315;
    
    if(control.generate_stress || control.stress_regime.empty()){
      prepare_stress();
    }
    else{
      stress_regime = control.stress_regime;
    }
    
  };

  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max) {
    std::cout << "Compute Environment: Time     "   << time << std::endl;
    const double lower_bound = 0.0;
    double upper_bound = height_max;

    auto f_canopy_openness = [&] (double height) -> double {return exp(-k_I * f_compute_competition(height));};
    environment_interpolator =
      environment_generator.construct(f_canopy_openness, lower_bound, upper_bound);
  }

  template <typename Function>
  void rescale_environment(Function f_compute_competition, double height_max) {
    std::vector<double> h = environment_interpolator.get_x();
    const double min = environment_interpolator.min(), // 0.0?
      height_max_old = environment_interpolator.max();

    auto f_canopy_openness = [&] (double height) -> double {return exp(-k_I * f_compute_competition(height));};
    util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
    h.back() = height_max; // Avoid round-off error.

    environment_interpolator.clear();
    for (auto hi : h) {
      environment_interpolator.add_point(hi, f_canopy_openness(hi));
    }
    environment_interpolator.initialise();
  }

  void set_fixed_environment(double value, double height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<double> y = {value, value, value};
    clear_environment();
    environment_interpolator.init(x, y);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double canopy_openness(double height) const {
    return get_environment_at_height(height);
  }
  
  
  double time_in_year() const{
    return (time-floor(time));
  }
  
  bool stressed() const {
    double yr = floor(time);
    if(time < stress_regime[yr]){
      return false;
    }
    else{
      return true;
    }
  }
  
  double getStress() const {
    double yr = floor(time);
    double yr_next = yr+1;
    double tcrit = stress_regime[yr];
    if(time < (yr + tcrit)/2){
      return 1/(1 + exp(-k_s * (time - yr)));
    }
    else if(time < (tcrit+yr_next)/2){
      return 1/(1 + exp(k_s * (time - tcrit)));
    }
    else{
      return 1/(1 + exp(-k_s * (time - yr_next)));
    }
  }
  
  double getStress_at_time(double t) const {
    double yr = floor(t);
    double yr_next = yr+1;
    double tcrit = stress_regime[yr];
    if(t < (yr + tcrit)/2){
      return 1/(1 + exp(-k_s * (t - yr)));
    }
    else if(t < (tcrit+yr_next)/2){
      return 1/(1 + exp(k_s * (t - tcrit)));
    }
    else{
      return 1/(1 + exp(-k_s * (t - yr_next)));
    }
  }
  
  void prepare_stress(){
    //erase stress as it exists
    stress_regime.erase(stress_regime.begin(),stress_regime.end());
    std::default_random_engine stress_regime_engine;

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<double> d(stress_mean,stress_sd);

    //hard coded 3000 as an upper limit that will probably not be reached
    for (int i = 0; i < 3000; i++) {
      stress_regime.push_back((d(stress_regime_engine) + i));
    }
  }

  void reset_stress_random(double new_mean=0.684931506849315, double new_sd=0.0684931506849315){
    stress_mean = new_mean;
    stress_sd = new_sd;
    prepare_stress();
  }

  void reset_stress(std::vector<double> new_stress_regime){
    stress_regime = new_stress_regime;
  }

  double stress_mean;
  double stress_sd;
  double k_s; // stress parameter
  std::vector<double> stress_regime;
};

inline Rcpp::NumericMatrix get_state(const ES20_Environment environment) {
  using namespace Rcpp;
  NumericMatrix xy = environment.environment_interpolator.r_get_xy();
  Rcpp::CharacterVector colnames =
    Rcpp::CharacterVector::create("height", "canopy_openness");
  xy.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);
  return xy;
}


}

#endif
