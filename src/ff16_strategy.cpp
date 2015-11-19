#include <plant/ff16_strategy.h>
#include <plant/uniroot.h>
#include <plant/qag.h>
#include <plant/environment.h>
#include <RcppCommon.h> // NA_REAL

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?

FF16_Strategy::FF16_Strategy() {
  // * Core traits - default values
  lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  rho       = 608.0;      // wood density [kg/m3]
  hmat      = 16.5958691; // Height at maturation [m]
  omega     = 3.8e-5;     // Seed mass [kg]

  // * Individual allometry
  // Canopy shape parameter (extra calculation here later)
  eta       = 12.0;
  // ratio leaf area to sapwood area
  theta     = 4669;
  // Height - leaf mass scaling
  a_l1        = 5.44;
  a_l2        = 0.306;
  // Root - leaf scaling
  a_r1        = 0.07;
  // Ratio of bark area : sapwood area
  a_b1         = 0.17;

  // * Production
  // Ratio of leaf dark respiration to leaf mass
  // =  [mol CO2 / kg(leaf) / yr]  |  (39.27)         |
  //    / [kg(leaf) / m2 ]         |    / (0.1978791) | lma
  // Hard coded in value of lma here so that this value doesn't change
  // if that trait changes above.
  r_l   = 39.27 / 0.1978791;
  // Root respiration per mass [mol CO2 / kg / yr]
  r_r   = 217.0;
  // Sapwood respiration per stem mass
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
  // Leaf turnover
  k_l=  0.4565855;
 // Bark turnover
  k_b    = 0.2;
  // Sapwood turnover
  k_s     = 0.2;
  // Root turnover
  k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  a_p1   = 151.177775377968;
  a_p2   = 0.204716166503633;

  // * Seed production
  // Accessory cost of reproduction, kg per seed
  a_f3  = 3.0 *  3.8e-5;

  // Maximum allocation to reproduction
  a_f1   = 1.0;
  // Size range across which individuals mature
  a_f2   = 50;

  // * Mortality parameters
  // Probability of survival during dispersal
  S_D   = 0.25;
  // Parameter for seedling survival
  a_d0    = 0.1;
  // Baseline for intrinsic mortality
  d_I    = 0.01;
 // Baseline rate for growth-related mortality
  a_dG1    = 5.5;
  // Risk coefficient for dry mass production (per area)
  a_dG2    = 20.0;

  // Will get computed properly by prepare_strategy
  height_0 = NA_REAL;
  eta_c    = NA_REAL;
}

// [eqn 2] area_leaf (inverse of [eqn 3])
double FF16_Strategy::area_leaf(double height) const {
  return pow(height / a_l1, 1.0 / a_l2);
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
double FF16_Strategy::mass_leaf(double area_leaf) const {
  return area_leaf * lma;
}

// [eqn 4] area and mass of sapwood
double FF16_Strategy::area_sapwood(double area_leaf) const {
  return area_leaf / theta;
}

double FF16_Strategy::mass_sapwood(double area_sapwood, double height) const {
  return area_sapwood * height * eta_c * rho;
}

// [eqn 5] area and mass of bark
double FF16_Strategy::area_bark(double area_leaf) const {
  return a_b1 * area_leaf / theta;
}

double FF16_Strategy::mass_bark(double area_bark, double height) const {
  return area_bark * height * eta_c * rho;
}

double FF16_Strategy::area_stem(double area_bark, double area_sapwood,
                            double area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

double FF16_Strategy::diameter_stem(double area_stem) const {
  return std::sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
double FF16_Strategy::mass_root(double area_leaf) const {
  return a_r1 * area_leaf;
}

// [eqn 8] Total mass
double FF16_Strategy::mass_live(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

double FF16_Strategy::mass_total(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_heartwood,
                            double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

double FF16_Strategy::mass_above_ground(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_root;
}

// one-shot update of the scm variables
void FF16_Strategy::scm_vars(const Environment& environment,
                              bool reuse_intervals,
                              Plant_internals& vars) {
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, vars.height, vars.area_leaf,
                           reuse_intervals);
  if (net_mass_production_dt_ > 0) {
    const double fraction_allocation_reproduction_ =
      fraction_allocation_reproduction(vars.height);
    const double darea_leaf_dmass_live_ =
      darea_leaf_dmass_live(vars.area_leaf);
    const double fraction_allocation_growth_ =
      fraction_allocation_growth(vars.height);
    const double area_leaf_dt =
      net_mass_production_dt_ * fraction_allocation_growth_ *
      darea_leaf_dmass_live_;
    vars.height_dt =
      dheight_darea_leaf(vars.area_leaf) * area_leaf_dt;
    vars.fecundity_dt =
      fecundity_dt(net_mass_production_dt_,
                   fraction_allocation_reproduction_);

    vars.area_heartwood_dt = area_heartwood_dt(vars.area_leaf);
    const double area_sapwood_ = area_sapwood(vars.area_leaf);
    const double mass_sapwood_ = mass_sapwood(area_sapwood_, vars.height);
    vars.mass_heartwood_dt = mass_heartwood_dt(mass_sapwood_);
  } else {
    vars.height_dt         = 0.0;
    vars.fecundity_dt      = 0.0;
    vars.area_heartwood_dt = 0.0;
    vars.mass_heartwood_dt = 0.0;
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.mortality_dt =
      mortality_dt(net_mass_production_dt_ / vars.area_leaf, vars.mortality);
}

// [eqn 12] Gross annual CO2 assimilation
//
// NOTE: In contrast with Daniel's implementation (but following
// Falster 2012), we do not normalise by a_y*a_bio here.
double FF16_Strategy::assimilation(const Environment& environment,
                                    double height,
                                    double area_leaf,
                                    bool reuse_intervals) {
  const bool over_distribution = control.plant_assimilation_over_distribution;
  const double x_min = 0.0, x_max = over_distribution ? 1.0 : height;

  double A = 0.0;

  std::function<double(double)> f;
  if (over_distribution) {
    f = [&] (double x) -> double {
      return compute_assimilation_p(x, height, environment);
    };
  } else {
    f = [&] (double x) -> double {
      return compute_assimilation_h(x, height, environment);
    };
  }

  if (control.plant_assimilation_adaptive && reuse_intervals) {
    A = control.integrator.integrate_with_last_intervals(f, x_min, x_max);
  } else {
    A = control.integrator.integrate(f, x_min, x_max);
  }

  return area_leaf * A;
}

// This is used in the calculation of assimilation by
// `compute_assimilation` above; it is the term within the integral in
// [eqn 12]; i.e., A_lf(A_0v, E(z,a)) * q(z,h(m_l))
// where `z` is height.
double FF16_Strategy::compute_assimilation_x(double x, double height,
                                     const Environment& environment) const {
  if (control.plant_assimilation_over_distribution) {
    return compute_assimilation_p(x, height, environment);
  } else {
    return compute_assimilation_h(x, height, environment);
  }
}

double FF16_Strategy::compute_assimilation_h(double z, double height,
                                     const Environment& environment) const {
  return assimilation_leaf(environment.canopy_openness(z)) * q(z, height);
}

double FF16_Strategy::compute_assimilation_p(double p, double height,
                                     const Environment& environment) const {
  return assimilation_leaf(environment.canopy_openness(Qp(p, height)));
}

// [Appendix S6] Per-leaf photosynthetic rate.
// Here, `x` is openness, ranging from 0 to 1.
double FF16_Strategy::assimilation_leaf(double x) const {
  return a_p1 * x / (x + a_p2);
}

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by a_y*a_bio.
double FF16_Strategy::respiration(double mass_leaf, double mass_sapwood,
                             double mass_bark, double mass_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root);
}

double FF16_Strategy::respiration_leaf(double mass) const {
  return r_l * mass;
}

double FF16_Strategy::respiration_bark(double mass) const {
  return r_b * mass;
}

double FF16_Strategy::respiration_sapwood(double mass) const {
  return r_s * mass;
}

double FF16_Strategy::respiration_root(double mass) const {
  return r_r * mass;
}

// [eqn 14] Total turnover
double FF16_Strategy::turnover(double mass_leaf, double mass_bark,
                          double mass_sapwood, double mass_root) const {
   return turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root);
}

double FF16_Strategy::turnover_leaf(double mass) const {
  return k_l * mass;
}

double FF16_Strategy::turnover_bark(double mass) const {
  return k_b * mass;
}

double FF16_Strategy::turnover_sapwood(double mass) const {
  return k_s * mass;
}

double FF16_Strategy::turnover_root(double mass) const {
  return k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
double FF16_Strategy::net_mass_production_dt_A(double assimilation, double respiration,
                                double turnover) const {
  return a_bio * a_y * (assimilation - respiration) - turnover;
}

// One shot calculation of net_mass_production_dt
// Used by germination_probability() and scm_vars().
double FF16_Strategy::net_mass_production_dt(const Environment& environment,
                                double height, double area_leaf_,
                                bool reuse_intervals) {
  const double mass_leaf_    = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);
  const double assimilation_ = assimilation(environment, height,
                                            area_leaf_, reuse_intervals);
  const double respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
    turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// [eqn 16] Fraction of production allocated to reproduction
double FF16_Strategy::fraction_allocation_reproduction(double height) const {
  return a_f1 / (1.0 + exp(a_f2 * (1.0 - height / hmat)));
}

// Fraction of production allocated to growth
double FF16_Strategy::fraction_allocation_growth(double height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
double FF16_Strategy::fecundity_dt(double net_mass_production_dt,
                               double fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (omega + a_f3);
}

double FF16_Strategy::darea_leaf_dmass_live(double area_leaf) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf));
}

// TODO: Ordering below here needs working on, probably as @dfalster
// does equation documentation?
double FF16_Strategy::dheight_darea_leaf(double area_leaf) const {
  return a_l1 * a_l2 * pow(area_leaf, a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
double FF16_Strategy::dmass_leaf_darea_leaf(double /* area_leaf */) const {
  return lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
double FF16_Strategy::dmass_sapwood_darea_leaf(double area_leaf) const {
  return rho * eta_c * a_l1 / (theta) * (a_l2 + 1.0) * pow(area_leaf, a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
double FF16_Strategy::dmass_bark_darea_leaf(double area_leaf) const {
  return a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
double FF16_Strategy::dmass_root_darea_leaf(double /* area_leaf */) const {
  return a_r1;
}

// Growth rate of basal diameter_stem per unit time
double FF16_Strategy::ddiameter_stem_darea_stem(double area_stem) const {
  return pow(M_PI / area_stem, 0.5);
}

// Growth rate of sapwood area at base per unit time
double FF16_Strategy::area_sapwood_dt(double area_leaf_dt) const {
  return area_leaf_dt / theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
double FF16_Strategy::area_heartwood_dt(double area_leaf) const {
  return k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
double FF16_Strategy::area_bark_dt(double area_leaf_dt) const {
  return a_b1 * area_leaf_dt / theta;
}

// Growth rate of stem basal area per unit time
double FF16_Strategy::area_stem_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
double FF16_Strategy::diameter_stem_dt(double area_stem, double area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
double FF16_Strategy::mass_root_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

double FF16_Strategy::mass_live_dt(double fraction_allocation_reproduction,
                               double net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

// TODO: Change top two to use mass_live_dt
double FF16_Strategy::mass_total_dt(double fraction_allocation_reproduction,
                                     double net_mass_production_dt,
                                     double mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
double FF16_Strategy::mass_above_ground_dt(double area_leaf,
                                       double fraction_allocation_reproduction,
                                       double net_mass_production_dt,
                                       double mass_heartwood_dt,
                                       double area_leaf_dt) const {
  const double mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

double FF16_Strategy::mass_heartwood_dt(double mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


double FF16_Strategy::mass_live_given_height(double height) const {
  double area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

double FF16_Strategy::height_given_mass_leaf(double mass_leaf) const {
  return a_l1 * pow(mass_leaf / lma, a_l2);
}

double FF16_Strategy::mortality_dt(double productivity_area,
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

double FF16_Strategy::mortality_growth_independent_dt() const {
  return d_I;
}

double FF16_Strategy::mortality_growth_dependent_dt(double productivity_area) const {
  return a_dG1 * exp(-a_dG2 * productivity_area);
}

// [eqn 20] Survival of seedlings during germination
double FF16_Strategy::germination_probability(const Environment& environment) {
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0);
  } else {
    return 0.0;
  }
}

double FF16_Strategy::area_leaf_above(double z, double height,
                                 double area_leaf) const {
  return area_leaf * Q(z, height);
}

// [eqn  9] Probability density of leaf area at height `z`
double FF16_Strategy::q(double z, double height) const {
  const double tmp = pow(z / height, eta);
  return 2 * eta * (1 - tmp) * tmp / z;
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height'
double FF16_Strategy::Q(double z, double height) const {
  if (z > height) {
    return 0.0;
  }
  const double tmp = 1.0-pow(z / height, eta);
  return tmp * tmp;
}

// (inverse of [eqn 10]; return the height above which fraction 'x' of
// the leaf mass would be found).
double FF16_Strategy::Qp(double x, double height) const { // x in [0,1], unchecked.
  return pow(1 - sqrt(x), (1/eta)) * height;
}

// The aim is to find a plant height that gives the correct seed mass.
double FF16_Strategy::height_seed(void) const {

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

void FF16_Strategy::prepare_strategy() {
  // Set up the integrator
  control.initialize();
  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2/(1 + eta) + 1/(1 + 2*eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
}

FF16_Strategy::ptr make_strategy_ptr(FF16_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16_Strategy>(s);
}

}
