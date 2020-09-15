// Comment this out for now: this will instead be a component like Canopy
// #include <plant/models/water_strategy.h>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
/*Water_Strategy::Water_Strategy() {
  // * Core traits - default values
  lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  rho       = 608.0;      // Wood density [kg/m3]
  hmat      = 16.5958691; // Height at maturation [m]
  omega     = 3.8e-5;     // Seed mass [kg]

  // * Individual allometry
  // Canopy shape parameter (extra calculation here later)
  eta       = 12.0; // [dimensionless]
  // Ratio sapwood area area to leaf area
  theta     = 1.0/4669; // [dimensionless]
  // Height - leaf mass scaling
  a_l1        = 5.44; // height with 1m2 leaf [m]
  a_l2        = 0.306; // dimensionless scaling of height with leaf area
  // Root mass per leaf area
  a_r1        = 0.07;  //[kg / m]
  // Ratio of bark area : sapwood area
  a_b1         = 0.17; // [dimensionless]

  // * Production
  // Ratio of leaf dark respiration to leaf mass [mol CO2 / yr  / kg]
  // =  [mol CO2 / m2 / yr]  |  (39.27 = 2100 * 0.00187)  | narea * photosynthesis_per_nitrogen
  //    / [kg(leaf) / m2 ]   |    / (0.1978791)           | lma
  // Hard coded in value of lma here so that this value doesn't change
  // if that trait changes above.
  r_l   = 39.27 / 0.1978791;
  // Root respiration per mass [mol CO2 / yr / kg]
  r_r   = 217.0;
  // Sapwood respiration per stem mass  [mol CO2 / yr / kg]
  // = respiration per volume [mol CO2 / m3 / yr]
  // /  wood density [kg/m3]
  r_s   = 4012.0 / 608.0;
  // Bark respiration per stem mass
  // assumed to be twice rate of sapwood
  // (NOTE that there is a re-parametrisation here relative to the paper
  // -- r_b is defined (new) as 2*r_s, whereas the paper assumes a
  // fixed multiplication by 2)
  r_b   = 2.0 * r_s;
  // Carbon conversion parameter
  a_y      = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  a_bio  = 2.45e-2;
  // Leaf turnover [/yr]
  k_l=  0.4565855;
 // Bark turnover [/yr]
  k_b    = 0.2;
  // Sapwood turnover [/yr]
  k_s     = 0.2;
  // Root turnover [/yr]
  k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  a_p1   = 151.177775377968; // [mol CO2 / yr / m2]
  a_p2   = 0.204716166503633; // [dimensionless]

  // * Seed production
  // Accessory cost of reproduction
  a_f3  = 3.0 *  3.8e-5; // [kg per seed]

  // Maximum allocation to reproduction
  a_f1   = 1.0; //[dimensionless]
  // Size range across which individuals mature
  a_f2   = 50; // [dimensionless]

// * Mortality parameters
  // Probability of survival during dispersal
  S_D   = 0.25; // [dimensionless]
  // Parameter for seedling survival
  a_d0    = 0.1; //[kg / yr / m2]
  // Baseline for intrinsic mortality
  d_I    = 0.01; // [ / yr]
 // Baseline rate for growth-related mortality
  a_dG1    = 5.5; // [ / yr]
  // Risk coefficient for dry mass production (per area)
  a_dG2    = 20.0;// [yr m2 / kg ]

  // Will get computed properly by prepare_strategy
  height_0 = NA_REAL;
  eta_c    = NA_REAL;

  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "Water";
}

void Water_Strategy::refresh_indices () {
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

// [eqn 2] area_leaf (inverse of [eqn 3])
double Water_Strategy::area_leaf(double height) const {
  return pow(height / a_l1, 1.0 / a_l2);
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
double Water_Strategy::mass_leaf(double area_leaf) const {
  return area_leaf * lma;
}

// [eqn 4] area and mass of sapwood
double Water_Strategy::area_sapwood(double area_leaf) const {
  return area_leaf * theta;
}

double Water_Strategy::mass_sapwood(double area_sapwood, double height) const {
  return area_sapwood * height * eta_c * rho;
}

// [eqn 5] area and mass of bark
double Water_Strategy::area_bark(double area_leaf) const {
  return a_b1 * area_leaf * theta;
}

double Water_Strategy::mass_bark(double area_bark, double height) const {
  return area_bark * height * eta_c * rho;
}

double Water_Strategy::area_stem(double area_bark, double area_sapwood,
                            double area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

double Water_Strategy::diameter_stem(double area_stem) const {
  return std::sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
double Water_Strategy::mass_root(double area_leaf) const {
  return a_r1 * area_leaf;
}

// [eqn 8] Total mass
double Water_Strategy::mass_live(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

double Water_Strategy::mass_total(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_heartwood,
                            double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

double Water_Strategy::mass_above_ground(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_root;
}

// for updating auxillary state
void Water_Strategy::update_dependent_aux(const int index, Internals& vars) {
  if (index == HEIGHT_INDEX) {
    double height = vars.state(HEIGHT_INDEX);
    vars.set_aux(aux_index.at("competition_effect"), area_leaf(height));
  }
}


// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void Water_Strategy::compute_rates(const Water_Environment& environment,
                              bool reuse_intervals,
                              Internals& vars) {

  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.aux(aux_index.at("competition_effect"));

  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_, reuse_intervals);
  // const double soil_water_use_dt_ =
    // soil_water_use_dt(environment, height, area_leaf_, reuse_intervals);

  // store the aux sate
  vars.set_aux(aux_index.at("net_mass_production_dt"), net_mass_production_dt_);
  //vars.set_aux(aux_index.at("soil_water_use_dt"), soil_water_use_dt_);

  if (net_mass_production_dt_ > 0) {
    
    const double fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
    const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_);
    const double fraction_allocation_growth_ = fraction_allocation_growth(height);
    const double area_leaf_dt = net_mass_production_dt_ * fraction_allocation_growth_ * darea_leaf_dmass_live_;
      
    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_) * area_leaf_dt);
    vars.set_rate(FECUNDITY_INDEX,
      fecundity_dt(net_mass_production_dt_, fraction_allocation_reproduction_));

    vars.set_rate(state_index.at("area_heartwood"), area_heartwood_dt(area_leaf_));
    const double area_sapwood_ = area_sapwood(area_leaf_);
    const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
    vars.set_rate(state_index.at("mass_heartwood"), mass_heartwood_dt(mass_sapwood_));

    if (collect_all_auxillary) {
      vars.set_aux(aux_index.at("area_sapwood"), area_sapwood_);
    }
  } else {
    vars.set_rate(HEIGHT_INDEX, 0.0);
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_rate(state_index.at("area_heartwood"), 0.0);
    vars.set_rate(state_index.at("mass_heartwood"), 0.0);
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(net_mass_production_dt_ / area_leaf_, vars.state(MORTALITY_INDEX)));
}


// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by a_y*a_bio.
double Water_Strategy::respiration(double mass_leaf, double mass_sapwood,
                             double mass_bark, double mass_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root);
}

double Water_Strategy::respiration_leaf(double mass) const {
  return r_l * mass;
}

double Water_Strategy::respiration_bark(double mass) const {
  return r_b * mass;
}

double Water_Strategy::respiration_sapwood(double mass) const {
  return r_s * mass;
}

double Water_Strategy::respiration_root(double mass) const {
  return r_r * mass;
}

// [eqn 14] Total turnover
double Water_Strategy::turnover(double mass_leaf, double mass_bark,
                          double mass_sapwood, double mass_root) const {
   return turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root);
}

double Water_Strategy::turnover_leaf(double mass) const {
  return k_l * mass;
}

double Water_Strategy::turnover_bark(double mass) const {
  return k_b * mass;
}

double Water_Strategy::turnover_sapwood(double mass) const {
  return k_s * mass;
}

double Water_Strategy::turnover_root(double mass) const {
  return k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
double Water_Strategy::net_mass_production_dt_A(double assimilation, double respiration,
                                double turnover) const {
  return a_bio * a_y * (assimilation - respiration) - turnover;
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
double Water_Strategy::net_mass_production_dt(const Water_Environment& environment,
                                double height, double area_leaf_,
                                bool reuse_intervals) {
  const double mass_leaf_    = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);
  const double assimilation_ = assimilator.assimilate(control, environment, height,
                                            area_leaf_, reuse_intervals);
  const double respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
    turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// [eqn 16] Fraction of production allocated to reproduction
double Water_Strategy::fraction_allocation_reproduction(double height) const {
  return a_f1 / (1.0 + exp(a_f2 * (1.0 - height / hmat)));
}

// Fraction of production allocated to growth
double Water_Strategy::fraction_allocation_growth(double height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
double Water_Strategy::fecundity_dt(double net_mass_production_dt,
                               double fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (omega + a_f3);
}

double Water_Strategy::darea_leaf_dmass_live(double area_leaf) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf));
}

// TODO: Ordering below here needs working on, probably as @dfalster
// does equation documentation?
double Water_Strategy::dheight_darea_leaf(double area_leaf) const {
  return a_l1 * a_l2 * pow(area_leaf, a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
double Water_Strategy::dmass_leaf_darea_leaf(double area_leaf) const {
  return lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
double Water_Strategy::dmass_sapwood_darea_leaf(double area_leaf) const {
  return rho * eta_c * a_l1 * theta * (a_l2 + 1.0) * pow(area_leaf, a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
double Water_Strategy::dmass_bark_darea_leaf(double area_leaf) const {
  return a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
double Water_Strategy::dmass_root_darea_leaf(double area_leaf) const {
  return a_r1;
}

// Growth rate of basal diameter_stem per unit time
double Water_Strategy::ddiameter_stem_darea_stem(double area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
double Water_Strategy::area_sapwood_dt(double area_leaf_dt) const {
  return area_leaf_dt * theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
double Water_Strategy::area_heartwood_dt(double area_leaf) const {
  return k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
double Water_Strategy::area_bark_dt(double area_leaf_dt) const {
  return a_b1 * area_leaf_dt * theta;
}

// Growth rate of stem basal area per unit time
double Water_Strategy::area_stem_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
double Water_Strategy::diameter_stem_dt(double area_stem, double area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
double Water_Strategy::mass_root_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

double Water_Strategy::mass_live_dt(double fraction_allocation_reproduction,
                               double net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

// TODO: Change top two to use mass_live_dt
double Water_Strategy::mass_total_dt(double fraction_allocation_reproduction,
                                     double net_mass_production_dt,
                                     double mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
double Water_Strategy::mass_above_ground_dt(double area_leaf,
                                       double fraction_allocation_reproduction,
                                       double net_mass_production_dt,
                                       double mass_heartwood_dt,
                                       double area_leaf_dt) const {
  const double mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

double Water_Strategy::mass_heartwood_dt(double mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


double Water_Strategy::mass_live_given_height(double height) const {
  double area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

double Water_Strategy::height_given_mass_leaf(double mass_leaf) const {
  return a_l1 * pow(mass_leaf / lma, a_l2);
}

double Water_Strategy::mortality_dt(double productivity_area,
                              double cumulative_mortality) const {

  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (R_FINITE(cumulative_mortality)) {
    return
      mortality_growth_independent_dt() +
      mortality_growth_dependent_dt(productivity_area);
 } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    return 0.0;
  }
}

double Water_Strategy::mortality_growth_independent_dt() const {
  return d_I;
}

double Water_Strategy::mortality_growth_dependent_dt(double productivity_area) const {
  return a_dG1 * exp(-a_dG2 * productivity_area);
}

// [eqn 20] Survival of seedlings during establishment
double Water_Strategy::establishment_probability(const Water_Environment& environment) {
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0);
  } else {
    return 0.0;
  }
}

double Water_Strategy::area_leaf_above(double z, double height) const {
  return area_leaf(height) * Q(z, height);
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height'
double Water_Strategy::Q(double z, double height) const {
  if (z > height) {
    return 0.0;
  }
  const double tmp = 1.0-pow(z / height, eta);
  return tmp * tmp;
}

// The aim is to find a plant height that gives the correct seed mass.
double Water_Strategy::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const double
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(omega);

  const double tol = control.plant_seed_tol;
  const size_t max_iterations = control.plant_seed_iterations;

  auto target = [&] (double x) mutable -> double {
    return mass_live_given_height(x) - omega;
  };

  return util::uniroot(target, h0, h1, tol, max_iterations);
}

void Water_Strategy::prepare_strategy() {
  // Set up the integrator
  control.initialize();
  assimilator.initialize(a_p1, a_p2, eta);
  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2/(1 + eta) + 1/(1 + 2*eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
}

Water_Strategy::ptr make_strategy_ptr(Water_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<Water_Strategy>(s);
}
*/
}
