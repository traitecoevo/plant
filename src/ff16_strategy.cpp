#include <plant/models/ff16_strategy.h>

namespace plant {

FF16_Strategy::FF16_Strategy() {
  collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16";
}

void FF16_Strategy::refresh_indices () {
    // Create and fill the name to state index maps
  state_index = std::map<std::string,int>();
  aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (size_t i = 0; i < state_names_vec.size(); i++) {
    state_index[state_names_vec[i]] = i;
  }
  for (size_t i = 0; i < aux_names_vec.size(); i++) {
    aux_index[aux_names_vec[i]] = i;
  }
}

// area_leaf() is defined inline in ff16_strategy.h (hot path).

// [eqn 1] mass_leaf (inverse of [eqn 2])
double FF16_Strategy::mass_leaf(double area_leaf) const {
  return area_leaf * lma;
}

// [eqn 4] area and mass of sapwood
double FF16_Strategy::area_sapwood(double area_leaf) const {
  return area_leaf * theta;
}

double FF16_Strategy::mass_sapwood(double area_sapwood, double height) const {
  return area_sapwood * height * eta_c * rho;
}

// [eqn 5] area and mass of bark
double FF16_Strategy::area_bark(double area_leaf) const {
  return a_b1 * area_leaf * theta;
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

// update_dependent_aux() is defined inline in ff16_strategy.h (hot path).

// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void FF16_Strategy::compute_rates(const FF16_Environment& environment,  Internals& vars) {

  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.aux(COMPETITION_EFFECT_AUX_INDEX);
  double height_inverse = vars.aux(HEIGHT_INVERSE_AUX_INDEX);

  // Reuse the sapwood intermediates the worker already computes (for
  // respiration/turnover) rather than recomputing them below; bit-identical.
  double area_sapwood_, mass_sapwood_;
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_, height_inverse,
                           area_sapwood_, mass_sapwood_);

  // store the aux sate
  vars.set_aux(NET_MASS_PRODUCTION_DT_AUX_INDEX, net_mass_production_dt_);

  if (net_mass_production_dt_ > 0) {

    const double fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
    // dheight_darea_leaf and the sapwood/bark terms in darea_leaf_dmass_live all
    // share pow(area_leaf, a_l2); evaluate this libm pow once and reuse it for
    // both rates rather than paying it twice per node per step (issue #361).
    const double area_leaf_pow_a_l2 = pow(area_leaf_, a_l2);
    const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_, area_leaf_pow_a_l2);
    const double fraction_allocation_growth_ = fraction_allocation_growth(height);
    const double area_leaf_dt = net_mass_production_dt_ * fraction_allocation_growth_ * darea_leaf_dmass_live_;

    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_, area_leaf_pow_a_l2) * area_leaf_dt);
    vars.set_rate(FECUNDITY_INDEX,
      fecundity_dt(net_mass_production_dt_, fraction_allocation_reproduction_));

    vars.set_rate(AREA_HEARTWOOD_INDEX, area_heartwood_dt(area_leaf_));
    vars.set_rate(MASS_HEARTWOOD_INDEX, mass_heartwood_dt(mass_sapwood_));

    if (collect_all_auxiliary) {
      vars.set_aux(AREA_SAPWOOD_AUX_INDEX, area_sapwood_);
    }
  } else {
    vars.set_rate(HEIGHT_INDEX, 0.0);
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_rate(AREA_HEARTWOOD_INDEX, 0.0);
    vars.set_rate(MASS_HEARTWOOD_INDEX, 0.0);
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(net_mass_production_dt_ / area_leaf_, vars.state(MORTALITY_INDEX)));
}

// [eqn 12] Gross annual CO2 assimilation -- deep-crown model.
// Integrate photosynthesis over crown depth: for a given height in the crown,
// take photosynthesis at that depth multiplied by the amount of leaf there.
double FF16_Strategy::assimilation_deep_crown(const FF16_Environment& environment,
                                              double height,
                                              double area_leaf,
                                              double height_inverse) {

  double A = 0.0;

  // Define an anonymous function to integrate.
  // Keep the lambda's own closure type (do not wrap in std::function) so the
  // templated QK::integrate inlines the integrand at each quadrature point
  // instead of making a type-erased indirect call.
  // Hoist the light-spline upper bound (canopy top) out of the integrand: it
  // is invariant across the quadrature, so fetch it once and pass it into the
  // capped get_environment_at_height() overload rather than re-reading
  // spline.max() at every quadrature point.
  const double canopy_top = environment.max_environment_height();
  auto f = [&](double z) -> double {
    return assimilation_leaf(environment.get_environment_at_height(z, canopy_top)) *
      canopy_shape.q(z * height_inverse, z);
  };

  // Integrate over crown depth using using Gauss-Kronrod quadrature.
  // The number of points used in the integration is determined by the control parameter
  // function_integration_rule. Rules defined in qk_rules.cpp
  A = function_integrator.integrate(f, 0.0, height);

  return area_leaf * A;
}

// [eqn 12] Gross annual CO2 assimilation -- mean-light model.
// Integrate the *light* over crown depth, weighted by the leaf-area density q
// (which integrates to one over the crown), to get the leaf-area-weighted mean
// light the crown experiences, then take a single photosynthesis evaluation of
// that mean. This sits between deep-crown (which integrates the concave
// photosynthetic rate itself) and crown-centre (a single point evaluation): it
// captures the mean light exactly but ignores the curvature of photosynthesis
// across the within-crown light distribution.
double FF16_Strategy::assimilation_average_light(const FF16_Environment& environment,
                                                 double height,
                                                 double area_leaf,
                                                 double height_inverse) {
  const double canopy_top = environment.max_environment_height();
  auto f = [&](double z) -> double {
    return environment.get_environment_at_height(z, canopy_top) *
      canopy_shape.q(z * height_inverse, z);
  };
  const double mean_light = function_integrator.integrate(f, 0.0, height);
  return area_leaf * assimilation_leaf(mean_light);
}

// [eqn 12] Gross annual CO2 assimilation -- crown-top model (crown-centre and PPA).
// Leaf area is treated as a thin layer at the crown centre, so a single
// evaluation of the light environment there replaces the crown-depth integral.
// The crown-centre and PPA models share this code: under crown-centre the environment
// returns the smooth light profile, under PPA it returns the stepped (layered)
// profile, so the only difference between them lives in the environment build.
// (height_inverse is unused but kept for a common dispatch signature.)
double FF16_Strategy::assimilation_crown_top(const FF16_Environment& environment,
                                             double height,
                                             double area_leaf,
                                             double /* height_inverse */) {
  const double E = environment.get_environment_at_height(height * eta_c);
  return area_leaf * assimilation_leaf(E);
}

// Photosynthetic rate per leaf area
// `x` is openness, ranging from 0 to 1.
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
// Used by establishment_probability() and compute_rates().
double FF16_Strategy::net_mass_production_dt(const FF16_Environment& environment,
                                double height, double area_leaf_) {
  return net_mass_production_dt(environment, height, area_leaf_,
                                1.0 / height);
}

double FF16_Strategy::net_mass_production_dt(const FF16_Environment& environment,
                                double height, double area_leaf_,
                                double height_inverse) {
  double area_sapwood_, mass_sapwood_;
  return net_mass_production_dt(environment, height, area_leaf_, height_inverse,
                                area_sapwood_, mass_sapwood_);
}

double FF16_Strategy::net_mass_production_dt(const FF16_Environment& environment,
                                double height, double area_leaf_,
                                double height_inverse,
                                double& area_sapwood_, double& mass_sapwood_) {
  const double mass_leaf_    = mass_leaf(area_leaf_);
  area_sapwood_ = area_sapwood(area_leaf_);
  mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);
  const double assimilation_ =
    assimilation(environment, height, area_leaf_, height_inverse);
  const double respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
    turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
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
  return darea_leaf_dmass_live(area_leaf, pow(area_leaf, a_l2));
}

double FF16_Strategy::darea_leaf_dmass_live(double area_leaf,
                                            double area_leaf_pow_a_l2) const {
  // dmass_bark_darea_leaf(area_leaf) == a_b1 * dmass_sapwood_darea_leaf(area_leaf),
  // so compute the shared pow(area_leaf, a_l2) term once rather than twice.
  const double dmass_sapwood_darea_leaf_ =
    dmass_sapwood_darea_leaf(area_leaf, area_leaf_pow_a_l2);
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf_
              + a_b1 * dmass_sapwood_darea_leaf_
              + dmass_root_darea_leaf(area_leaf));
}

double FF16_Strategy::dheight_darea_leaf(double area_leaf) const {
  return a_l1 * a_l2 * pow(area_leaf, a_l2 - 1);
}

double FF16_Strategy::dheight_darea_leaf(double area_leaf,
                                         double area_leaf_pow_a_l2) const {
  // pow(area_leaf, a_l2 - 1) == pow(area_leaf, a_l2) / area_leaf, so reuse the
  // shared power instead of a second libm pow (a reciprocal-multiply reorder,
  // so not bit-identical to the single-argument form).
  return a_l1 * a_l2 * area_leaf_pow_a_l2 / area_leaf;
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
double FF16_Strategy::dmass_leaf_darea_leaf(double /* area_leaf */) const {
  return lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
double FF16_Strategy::dmass_sapwood_darea_leaf(double area_leaf) const {
  return dmass_sapwood_darea_leaf(area_leaf, pow(area_leaf, a_l2));
}

double FF16_Strategy::dmass_sapwood_darea_leaf(double /* area_leaf */,
                                               double area_leaf_pow_a_l2) const {
  return rho * eta_c * a_l1 * theta * (a_l2 + 1.0) * area_leaf_pow_a_l2;
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
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
double FF16_Strategy::area_sapwood_dt(double area_leaf_dt) const {
  return area_leaf_dt * theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
double FF16_Strategy::area_heartwood_dt(double area_leaf) const {
  return k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
double FF16_Strategy::area_bark_dt(double area_leaf_dt) const {
  return a_b1 * area_leaf_dt * theta;
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

double FF16_Strategy::mass_total_dt(double fraction_allocation_reproduction,
                                     double net_mass_production_dt,
                                     double mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO(#480): Do we not track root mass change?
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
  if (util::is_finite(cumulative_mortality)) {
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

// [eqn 20] Survival of seedlings during establishment
double FF16_Strategy::establishment_probability(const FF16_Environment& environment) {
  
  double decay_over_time = exp(-recruitment_decay * environment.time);
  
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0,
                           height_0_inverse);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

// compute_competition() overloads and compute_competition_by_ratio() are
// defined inline in ff16_strategy.h (per-node hot path).

// (inverse of [eqn 10]; return the height above which fraction 'x' of
// the leaf mass would be found).
double FF16_Strategy::Qp(double x, double height) const { // x in [0,1], unchecked.
  return canopy_shape.Qp(x, height);
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

  const double tol = control.offspring_production_tol;
  const size_t max_iterations = control.offspring_production_iterations;

  auto target = [&] (double x) mutable -> double {
    return mass_live_given_height(x) - omega;
  };

  return util::uniroot(target, h0, h1, tol, max_iterations);
}

void FF16_Strategy::prepare_strategy() {

  // Set up the function_integrator
  function_integrator = quadrature::QK(
      // Gauss-Kronrod quadrature integeration rule (see qkrules)
      control.function_integration_rule);

  // Resolve the crown shading model once (string -> enum), then bind both hot
  // paths to it: canopy_shape handles competition (leaf_area_above), and
  // assimilation_fn selects the matching assimilation implementation. After
  // this, neither path compares the model string per call.
  const ShadingModel shading_model =
    shading_model_from_string(control.shading_model, ShadingModel::DeepCrown);
  // canopy_shape also selects the competition contribution: smooth Q for every
  // model except flat-top-box, which casts a step (see CanopyShape).
  canopy_shape.initialise(eta, shading_model);
  switch (shading_model) {
  case ShadingModel::DeepCrown:
    assimilation_fn = &FF16_Strategy::assimilation_deep_crown;
    break;
  case ShadingModel::MeanLight:
    assimilation_fn = &FF16_Strategy::assimilation_average_light;
    break;
  case ShadingModel::CrownCentre:
  case ShadingModel::FlatTopBox:
  case ShadingModel::FlatTopSoftBox:
  case ShadingModel::PPA:
    // All evaluate assimilation at the crown centre. They differ in the light
    // profile they read: crown-centre from the smooth profile, flat-top-box / -soft-
    // box from a profile built with (hard / smoothed) box competition, PPA from
    // a stepped profile.
    assimilation_fn = &FF16_Strategy::assimilation_crown_top;
    break;
  }

  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2/(1 + eta) + 1/(1 + 2*eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  height_0_inverse = 1.0 / height_0;
  area_leaf_0 = area_leaf(height_0);

  if (is_variable_birth_rate) {
    extrinsic_drivers.set_variable("birth_rate", birth_rate_x, birth_rate_y);
  } else {
    extrinsic_drivers.set_constant("birth_rate", birth_rate_y[0]);
  }
}

FF16_Strategy::ptr make_strategy_ptr(FF16_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16_Strategy>(s);
}
}
