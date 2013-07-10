#include "util.h"   // check_bounds_r

#include "parameters.h"

namespace model {

Parameters::Parameters() {
  reset();
}

void Parameters::reset() {
  // Light extinction coefficient
  c_ext = 0.5;
  // Patch area (m^2?)
  patch_area = 10.0;
  // Survival during dispersal
  Pi_0 = 0.25;
  // Number of patches
  n_patches = 10;
  // Mean disturbance interval in years
  mean_disturbance_interval = 30.0;

  // Then set the values for the lookup table, based on these (this is
  // basically the inverse of set_parameters_post_hook())
  _n_patches = static_cast<double>(n_patches);

  strategies.clear();
}

size_t Parameters::size() const {
  return strategies.size();
}

void Parameters::add_strategy(Strategy s) {
  s.set_control(control);
  strategies.push_back(s);
}

const Control& Parameters::get_control() const {
  return control;
}

void Parameters::set_control(Control x) {
  control = x;
  push_control_to_strategies();
}

void Parameters::set_control_parameters(Rcpp::List x) {
  control.set_parameters(x);
  push_control_to_strategies();
}

Strategy Parameters::r_at(size_t idx) {
  return strategies.at(util::check_bounds_r(idx, size()));
}

Rcpp::List Parameters::r_get_strategies() {
  Rcpp::List ret;
  for (strategy_iterator it = strategies.begin();
       it != strategies.end(); ++it)
    ret.push_back(Rcpp::wrap(*it));
  return ret;
}

void Parameters::do_build_lookup() {
  lookup_table["c_ext"]      = &c_ext;
  lookup_table["patch_area"] = &patch_area;
  lookup_table["Pi_0"]       = &Pi_0;
  lookup_table["n_patches"]  = &_n_patches;
  lookup_table["mean_disturbance_interval"] =
    &mean_disturbance_interval;
}

void Parameters::set_parameters_post_hook() {
  n_patches = static_cast<size_t>(_n_patches);
}

void Parameters::push_control_to_strategies() {
  for (strategy_iterator it = strategies.begin();
       it != strategies.end(); ++it)
    it->set_control(control);
}

}
