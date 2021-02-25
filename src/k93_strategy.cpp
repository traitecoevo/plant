// Built from  src/ff16r_strategy.cpp on Wed Aug 12 15:46:38 2020 using the scaffolder, from the strategy:  FF16r
// Built from  src/ff16_strategy.cpp on Fri Jul  3 08:14:35 2020 using the scaffolder, from the strategy:  FF16
#include <plant/uniroot.h>
#include <plant/qag.h>
#include <plant/models/assimilation.h>
#include <plant/models/ff16_environment.h>
#include <plant/models/k93_strategy.h>
#include <RcppCommon.h> // NA_REAL

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
    double height = vars.state(HEIGHT_INDEX);
    vars.set_aux(aux_index.at("competition_effect"),
                 compute_competition(0.0, height));
  }
}

// Smoothing function for competition effect
double K93_Strategy::Q(double z, double size) const {
  if (z > size) {
    return 0.0;
  }
  const double tmp = 1.0 - pow(z / size, eta);

  return tmp * tmp;
}

double K93_Strategy::compute_competition(double z, double size) const {

  // Competition only felt if plant bigger than target size z
  return size_to_basal_area(size) * Q(z, size);
 };

double K93_Strategy::establishment_probability(const K93_Environment& environment){
  //TODO: may want to make this dependent on achieving positive growth rate
  return 1.0;
}

double K93_Strategy::net_mass_production_dt(const K93_Environment& environment,
                                            double height, double area_leaf_,
                                            bool reuse_intervals) {
  // TODO: there was no return value here - added 0.0
  return 1.0;
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

  double height = vars.state(HEIGHT_INDEX);

  // suppression integral mapped [0, 1] using adaptive spline
  // back transform to basal area and add suppression from self
  double competition = environment.get_environment_at_height(height);
  double basal_area = size_to_basal_area(height);
  const double k_I = get_k_I(environment);

  double cumulative_basal_area = -log(competition) / k_I;

  if (!util::is_finite(cumulative_basal_area)) {
    util::stop("Environmental interpolator has gone out of bounds, try lowering the extinction coefficient k_I");
  }

  vars.set_rate(HEIGHT_INDEX,
    size_dt(height, cumulative_basal_area));

  vars.set_rate(FECUNDITY_INDEX,
    fecundity_dt(height, cumulative_basal_area));

  vars.set_rate(MORTALITY_INDEX,
    mortality_dt(cumulative_basal_area, vars.state(MORTALITY_INDEX)));
}

double K93_Strategy::size_to_basal_area(double size) const {
  return M_PI / 4 * pow(size, 2);
}

// [eqn 10] Growth
double K93_Strategy::size_dt(double size,
                             double cumulative_basal_area) const {

  double growth = size * (b_0 - b_1 * log(size) - b_2 * cumulative_basal_area);

  if(growth < 0.0) {
    growth = 0.0;
  }

  return growth;
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
    return (mu > 0)? mu:0.0;
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
