// Built from  src/ff16_strategy.cpp on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16
#include <plant/uniroot.h>
#include <plant/qag.h>
#include <plant/models/k93_strategy.h>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
K93_Strategy::K93_Strategy() {
   // * Empirical parameters - Table 1.
   height_0 = 2.0;
   b_0 = 0.059;    // Growth intercept year-1
   b_1 = 0.012;    // Growth asymptote year-1.(ln cm)-1
   b_2 = 0.00041;  // Growth suppression rate m2.cm-2.year-1
   c_0 = 0.008;    // Mortality intercept year-1
   c_1 = 0.00044;  // Mortality suppression rate m2.cm-2.year-1
   d_0 = 0.00073;  // Recruitment rate (cm2.year-1)
   d_1 = 0.044;    // Recruitment suppression rate (m2.cm-2)

  // build the string state/aux name to index map
  refresh_indices();
  name = "K93";
}

// Signatures fixed in plant.h
void K93_Strategy::update_dependent_aux(const int index, Internals& vars) {
  if (index == HEIGHT_INDEX) {
    double size = vars.state(HEIGHT_INDEX);
    vars.set_aux(aux_index.at("competition_effect"),
                 compute_competition(0.0, size));
  }
}

double K93_Strategy::establishment_probability(const K93_Environment& environment){
    return 1.0;
}

double K93_Strategy::net_mass_production_dt(const K93_Environment& environment,
                                            double height, double area_leaf_,
                                            bool reuse_intervals) {
}

void K93_Strategy::refresh_indices () {
  // Create and fill the name to state index maps
  state_index = std::map<std::string,int>();
  aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (int i = 0; i < state_names_vec.size(); i++) {
    state_index[state_names_vec[i]] = i;
  }
  for (int i = 0; i < aux_names_vec.size(); i++) {
    aux_index[aux_names_vec[i]] = i;
  }
}

// i.e. setting rates of ode vars from the state and updating aux vars
void K93_Strategy::compute_rates(const K93_Environment& environment,
                              bool reuse_intervals,
                              Internals& vars) {

  double size = vars.state(HEIGHT_INDEX);

  // suppression integral mapped [0, 1] using adaptive spline
  // back transform to basal area and add supression from self
  double competition = environment.get_environment_at_height(size);
  double basal_area = size_to_basal_area(size);

  double cumulative_basal_area = -log(competition) / environment.k_I + basal_area;

  vars.set_rate(HEIGHT_INDEX,
    size_dt(size, cumulative_basal_area));

  vars.set_rate(FECUNDITY_INDEX,
    fecundity_dt(size, cumulative_basal_area));

  vars.set_rate(MORTALITY_INDEX,
    mortality_dt(cumulative_basal_area, vars.state(MORTALITY_INDEX)));
}

double K93_Strategy::size_to_basal_area(double size) const {
  return M_PI / 4 * pow(size, 2);
}

// [eqn 10] Growth
double K93_Strategy::size_dt(double size,
                             double cumulative_basal_area) const {
  return pow(size, b_0 - b_1 * log(size) - b_2 * cumulative_basal_area);
}

// [eqn 12] Reproduction
double K93_Strategy::fecundity_dt(double size,
                                  double cumulative_basal_area) const {
  double basal_area = size_to_basal_area(size);
  return d_0 * basal_area * exp(-d_1 * cumulative_basal_area);
}

// [eqn 11] Mortality
double K93_Strategy::mortality_dt(double cumulative_basal_area,
                                  double cumulative_mortality) const {
  // If mortality probability is 1 (latency = Inf) then the rate
  // calculations break.  Setting them to zero gives the correct
  // behaviour.
  if (R_FINITE(cumulative_mortality)) {
    double mu = -c_0 + c_1 * cumulative_basal_area;
    return c_0 > mu + c_0 ? 0.0 : mu;
 } else {
    return 0.0;
  }
}

void K93_Strategy::prepare_strategy() {
  // Set up the integrator
  control.initialize();
}

K93_Strategy::ptr make_strategy_ptr(K93_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<K93_Strategy>(s);
}
}
