#include "plant.h"

#include "find_root.h"

namespace model {

Plant::Plant(Strategy s)
  : strategy(s) {
  set_mass_leaf(strategy->mass_leaf_0);
}

Plant::Plant(Strategy *s)
  : strategy(s) {
  set_mass_leaf(strategy->mass_leaf_0);
}

// TODO: Value of 0.0 not ideal, but allows more easy comparison.
// A value of NA_REAL would be better, but this requires that '=='
// won't work sensibly on freshly-produced plants.
Plant::internals::internals() 
  : mass_leaf(NA_REAL),
    leaf_area(NA_REAL),
    height(NA_REAL),
    mass_sapwood(NA_REAL),
    mass_bark(NA_REAL),
    mass_heartwood(NA_REAL),
    mass_root(NA_REAL),
    mass_total(NA_REAL),
    // TODO: Not sure about these being zero
    assimilation(0.0),
    respiration(0.0),
    turnover(0.0),
    net_production(0.0),
    reproduction_fraction(0.0),
    fecundity_rate(0.0),
    leaf_fraction(0.0),
    mass_leaf_growth_rate(0.0),
    height_growth_rate(0.0),
    mortality_rate(0.0),
    // But these should be zero
    mortality(0.0),
    fecundity(0.0) {}

// NOTE: This is essentially only used for testing.
// NOTE: The semantics around comparing Strategy are ill-defined for
// the standalone case.
bool Plant::operator==(const Plant &rhs) {
  if ( strategy.standalone() && strategy.get() == rhs.strategy.get() )
    Rf_warning("This is going to end badly.");
  return strategy == rhs.strategy && vars == rhs.vars;
}

bool Plant::internals::operator==(const Plant::internals &rhs) {
  const bool ret = true
    // Size members...
    && mass_leaf      == rhs.mass_leaf
    && leaf_area      == rhs.leaf_area
    && height         == rhs.height
    && mass_sapwood   == rhs.mass_sapwood
    && mass_bark      == rhs.mass_bark
    && mass_heartwood == rhs.mass_heartwood
    && mass_root      == rhs.mass_root
    && mass_total     == rhs.mass_total
    // // ...physiological members...
    && assimilation   == rhs.assimilation
    && respiration    == rhs.respiration
    && turnover       == rhs.turnover
    && net_production == rhs.net_production
    && reproduction_fraction == rhs.reproduction_fraction
    && fecundity_rate == rhs.fecundity_rate
    && leaf_fraction  == rhs.leaf_fraction
    && mass_leaf_growth_rate == rhs.mass_leaf_growth_rate
    && height_growth_rate == rhs.height_growth_rate
    && mortality_rate == rhs.mortality_rate
    // // ...and "other" members...
    && fecundity      == rhs.fecundity
    && mortality      == rhs.mortality
    ;  

  return ret;
}

// [eqn 1-8] Update size variables given an input leaf mass
// 
// NOTE: Only recomputes size variables if the mass is actually
// different.  This is totally safe if nothing else sets either
// mass_leaf or any size variable except this method.  This will help
// save quite a bit of calculation time and book-keeping down the
// track.
void Plant::set_mass_leaf(double mass_leaf_) {
  if ( mass_leaf_ <= 0.0 )
    Rf_error("mass_leaf must be positive (given %2.5f)", mass_leaf_);
  if ( mass_leaf_ != mass_leaf() )
    compute_vars_size(mass_leaf_);
}

// If we can set it, we can get it
double Plant::mass_leaf() const {
  return vars.mass_leaf;
}

double Plant::height() const {
  return vars.height;
}

void Plant::set_height(double height_) {
  set_mass_leaf(mass_leaf_given_height(height_));
}

double Plant::height_rate() const {
  return vars.height_growth_rate;
}

// * Competitive environment
// [      ] Leaf area (not fraction) above height `z`
double Plant::leaf_area_above(double z) const {
  return vars.leaf_area * Q(z);
}

// * Mass production
// [eqn 12-19,21] Update physiological variables given the current
// light environment (and given the current set of size variables).
void Plant::compute_vars_phys(const Environment& environment) {
  // [eqn 12] Gross annual CO2 assimilation
  vars.assimilation = compute_assimilation(environment);

  // [eqn 13] Total maintenance respiration
  vars.respiration = compute_respiration();

  // [eqn 14] Total turnover 
  vars.turnover = compute_turnover();

  // [eqn 15] Net production
  // 
  // NOTE: Translation of variable names from the EBT.  our
  // `net_primary_production` is EBT's N, our `net_production` is
  // EBT's P.
  const double net_primary_production = 
    strategy->c_bio * strategy->Y * (vars.assimilation - vars.respiration);
  vars.net_production = net_primary_production - vars.turnover;

  if ( vars.net_production > 0 ) {
    // [eqn 16] - Fraction of whole plant growth that is leaf
    vars.reproduction_fraction = compute_reproduction_fraction();

    // [eqn 17] - Rate of offspring production
    // 
    // NOTE: In EBT, was multiplied by Pi_0 (survival during
    // dispersal), but we do not here.
    vars.fecundity_rate = vars.net_production *
      vars.reproduction_fraction / (strategy->c_acc * strategy->s);

    // [eqn 18] - Fraction of mass growth in leaves
    vars.leaf_fraction = compute_leaf_fraction();

    // TODO: move into compute_growth_rate?

    // [eqn 19] - Growth rate in leaf mass
    vars.mass_leaf_growth_rate = vars.net_production *
      (1 - vars.reproduction_fraction) * vars.leaf_fraction;

    vars.height_growth_rate =
      strategy->a1 * strategy->B1 *
      pow(vars.leaf_area, strategy->B1 - 1) * vars.mass_leaf_growth_rate /
      strategy->lma;
  } else {
    vars.reproduction_fraction = 0.0;
    vars.fecundity_rate        = 0.0;
    vars.leaf_fraction         = 0.0;
    vars.mass_leaf_growth_rate = 0.0;
    vars.height_growth_rate    = 0.0;
  }

  // [eqn 21] - Instantaneous mortality rate
  //
  // Composed of a wood density effect (term involving c_d0) and a
  // growth effect (term involving c_d2)
  vars.mortality_rate = 
    strategy->c_d0 * exp(-strategy->c_d1 * strategy->rho) +
    strategy->c_d2 * exp(-strategy->c_d3 * 
			 vars.net_production / vars.leaf_area);
}

// * Births and deaths
int Plant::offspring() {
  double born = 0;
  if ( vars.fecundity > 1 )
    vars.fecundity = modf(vars.fecundity, &born);
  return (int)born;
}

bool Plant::died() {
  const int did_die = unif_rand() < mortality_probability();
  set_mortality(0.0);
  return did_die;
}

// [eqn 20] Survival of seedlings during germination
//
// NOTE: This does not check/enforce that mass_leaf is set to the seed
// mass leaf (so this is actually the germination probability of a
// plant that happens to be the current size).  This might be
// something to change.
double Plant::germination_probability(const Environment& environment) {
  compute_vars_phys(environment);
  if ( vars.net_production > 0 ) {
    const double tmp = vars.leaf_area * strategy->c_s0 / vars.net_production;
    return 1 / (tmp * tmp + 1.0);
  } else {
    return 0.0;
  }
}

// * Access the Control parameter
// NOTE: This might change!
const Control& Plant::control() const {
  return strategy->control;
}

// * ODE interface
size_t Plant::ode_size() const { 
  return ode_dimension; 
}

ode::iter_const Plant::ode_values_set(ode::iter_const it) {
  set_mass_leaf(*it++);
  set_mortality(*it++);
  set_fecundity(*it++);
  return it;
}

ode::iter Plant::ode_values(ode::iter it) const {
  *it++ = mass_leaf();
  *it++ = mortality();
  *it++ = fecundity();
  return it;
}

ode::iter Plant::ode_rates(ode::iter it) const {
  *it++ = mass_leaf_rate();
  *it++ = mortality_rate();
  *it++ = fecundity_rate();
  return it;
}

// * Protected methods

// These methods protected partly because they require that
// compute_vars_phys to be run before these will give sensible
// answers.

double Plant::mass_leaf_rate() const {
  return vars.mass_leaf_growth_rate;
}

double Plant::mortality() const {
  return vars.mortality;
}
void Plant::set_mortality(double x) {
  vars.mortality = x;
}
double Plant::mortality_rate() const {
  return vars.mortality_rate;
}

double Plant::fecundity() const {
  return vars.fecundity;
}
void Plant::set_fecundity(double x) {
  vars.fecundity = x;
}
double Plant::fecundity_rate() const {
  return vars.fecundity_rate;
}

// This one is a bit different, as it converts from the mean of the
// poisson process (on [0,Inf)) to a probability (on [0,1]).
double Plant::mortality_probability() const {
  return 1 - exp(-mortality());
}
double Plant::survival_probability() const {
  return exp(-mortality());
}

// * Private methods

// * Individual size
// [eqn 1-8] Update size variables to a new leaf mass.
void Plant::compute_vars_size(double mass_leaf_) {
  // [eqn 1] Leaf mass
  vars.mass_leaf = mass_leaf_;
  // [eqn 2] Leaf area
  vars.leaf_area = vars.mass_leaf / strategy->lma;
  // [eqn 3] Height
  vars.height = strategy->a1*pow(vars.leaf_area, strategy->B1);
  // [eqn 4] Mass of sapwood
  vars.mass_sapwood =   strategy->rho / strategy->theta *
    strategy->a1 * strategy->eta_c * pow(vars.leaf_area, 1 + strategy->B1);
  // [eqn 5] Mass of bark
  vars.mass_bark = strategy->b * vars.mass_sapwood;
  // [eqn 6] Mass of heartwood
  vars.mass_heartwood = strategy->rho * strategy->eta_c * strategy->a2 *
    pow(vars.leaf_area, strategy->B2);
  // [eqn 7] Mass of (fine) roots
  vars.mass_root = strategy->a3 * vars.leaf_area;
  // [eqn 8] Total mass
  vars.mass_total =
    vars.mass_leaf + vars.mass_sapwood + vars.mass_bark + 
    vars.mass_heartwood + vars.mass_root;
}


// [eqn  9] Probability density of leaf area at height `z`
double Plant::q(double z) const {
  const double eta = strategy->eta;
  const double tmp = pow(z / vars.height, eta);
  return 2 * eta * (1 - tmp) * tmp / z;
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height'
double Plant::Q(double z) const {
  if ( z > vars.height )
    return 0.0;
  const double tmp = 1.0-pow(z / vars.height, strategy->eta);
  return tmp * tmp;
}

// (inverse of [eqn 10]; return the height above which fraction 'x' of
// the leaf mass would be found).
double Plant::Qp(double x) const { // x in [0,1], unchecked.
  return pow(1 - sqrt(x), (1/strategy->eta)) * vars.height;
}

// [eqn 12] Gross annual CO2 assimilation
// 
// NOTE: In contrast with EBT, we do not normalise by Y*c_bio.
// 
// TODO: This version is completely naive.  Better to provide an
// already working integrator.
double Plant::compute_assimilation(const Environment& environment) const {
  FunctorBind1<Plant, const Environment&,
	       &Plant::compute_assimilation_x> fun(this, environment);
  const double tol = control().plant_assimilation_tol;
  const int max_iterations = control().plant_assimilation_iterations;
  util::Integrator integrator(tol, tol, max_iterations);
  const double x_max = control().plant_assimilation_over_distribution ?
    1 : vars.height;
  return vars.leaf_area * integrator.integrate(&fun, 0.0, x_max);
}

// This is used in the calculation of assimilation by
// `compute_assimilation` above; it is the term within the integral in
// [eqn 12]; i.e., A_lf(A_0v, E(z,a)) * q(z,h(m_l))
// where `z` is height.
double Plant::compute_assimilation_x(double x, 
				     const Environment& environment) const {
  if ( control().plant_assimilation_over_distribution )
    return assimilation_leaf(environment.canopy_openness(Qp(x)));
  else
    return assimilation_leaf(environment.canopy_openness(x)) * q(x);
}

// [Appendix S6] Per-leaf photosynthetic rate.
// Here, `x` is openness, ranging from 0 to 1.
double Plant::assimilation_leaf(double x) const {
  return strategy->c_p1 * x / (x + strategy->c_p2);
}

// [eqn 13] Total maintenance respiration
// 
// (NOTE that there is a reparametrisation here relative to the paper
// -- c_Rb is defined (new) as 2*c_Rs, wheras the paper assumes a
// fixed multiplication by 2)
//
// NOTE: In contrast with EBT, we do not normalise by Y*c_bio.
double Plant::compute_respiration() const {
  return
    strategy->c_Rl * vars.leaf_area * strategy->n_area +
    strategy->c_Rs * vars.mass_sapwood / strategy->rho +
    strategy->c_Rb * vars.mass_bark    / strategy->rho +
    strategy->c_Rr * vars.mass_root;
}

// [eqn 14] Total turnover
// 
// (NOTE: `k_l` is (a_4*\phi)^{b_4} in [eqn 14], and is computed by
// `prepare_strategy`).
double Plant::compute_turnover() const {
  return
    vars.mass_leaf * strategy->k_l  +
    vars.mass_bark * strategy->k_b  +
    vars.mass_root * strategy->k_r;
}

// [eqn 16] Fraction of whole plan growth that is leaf
double Plant::compute_reproduction_fraction() const {
  return strategy->c_r1 / (1.0 + exp(strategy->c_r2 *
				     (1.0 - vars.height/strategy->hmat)));
}

// [eqn 18] Fraction of mass growth that is leaves
// 
// TODO: This could do with documenting properly and tidying.
//
// NOTE: The EBT version actually computed 1/leaf_fraction (modifying
// growth rate calculation accordingly).  Possibly more stable?
double Plant::compute_leaf_fraction() const {
  const Strategy *s = strategy.get(); // for brevity.
  return 1.0/(1.0 + s->a3/s->lma +
	      (s->rho / s->theta * s->a1 * s->eta_c * (1.0 +s->b) *
	       (1.0+s->B1) * pow(vars.leaf_area, s->B1) / s->lma +
	       s->rho * s->a2 * s->eta_c * s->B2 *
	       pow(vars.leaf_area, s->B2-1) / s->lma));
}

// NOTE: static method
double Plant::mass_leaf_seed(Strategy *s) {
  Plant p(s);
  const double mass_seed = s->s;

  // Functor computing total mass of a seed with mass m, minus the
  // target seed mass.
  util::FunctorRoot<Plant, &Plant::compute_mass_total> fun(&p, mass_seed);

  const double tol = p.control().plant_seed_tol;
  const int max_iterations = p.control().plant_seed_iterations;
  util::RootFinder root(tol, tol, max_iterations);

  return root.root(&fun, DBL_MIN, mass_seed);
}

// NOTE: This is used only by mass_leaf_seed.
double Plant::compute_mass_total(double x) {
  set_mass_leaf(x);
  return vars.mass_total;
}

// NOTE: static method.
void Plant::prepare_strategy(Strategy *s) {
  s->eta_c = 1 - 2/(1 + s->eta) + 1/(1 + 2*s->eta);
  s->k_l = s->a4 * pow(s->lma, -s->B4);
  s->mass_leaf_0 = mass_leaf_seed(s);
}

// * R interface
Strategy Plant::r_get_strategy() const {
  return *strategy.get();
}

Rcpp::NumericVector Plant::r_get_vars_size() const {
  using namespace Rcpp;
  return NumericVector::create(_["mass_leaf"]=vars.mass_leaf,
			       _["mass_sapwood"]=vars.mass_sapwood,
			       _["mass_bark"]=vars.mass_bark,
			       _["mass_heartwood"]=vars.mass_heartwood,
			       _["mass_root"]=vars.mass_root,
			       _["mass_total"]=vars.mass_total,
			       _["height"]=vars.height,
			       _["leaf_area"]=vars.leaf_area);
}

Rcpp::NumericVector Plant::r_get_vars_phys() const {
  using namespace Rcpp;
  return NumericVector::create(_["assimilation"]=vars.assimilation,
			       _["respiration"]=vars.respiration,
			       _["turnover"]=vars.turnover,
			       _["net_production"]=vars.net_production,
			       _["reproduction_fraction"]=
			       vars.reproduction_fraction,
			       _["fecundity_rate"]=vars.fecundity_rate,
			       _["leaf_fraction"]=vars.leaf_fraction,
			       _["mass_leaf_growth_rate"]=
			       vars.mass_leaf_growth_rate,
			       _["height_growth_rate"]=
			       vars.height_growth_rate,
			       _["mortality_rate"]=vars.mortality_rate);
}

bool Plant::r_died() {
  Rcpp::RNGScope scope;
  return died();
}

double Plant::height_given_mass_leaf(double mass_leaf_) const {
  return strategy->a1 * pow(mass_leaf_ / strategy->lma, strategy->B1);
}

double Plant::mass_leaf_given_height(double height_) const {
  return strategy->lma * pow(height_ / strategy->a1, 1/strategy->B1);
}

namespace test {

bool test_plant(Strategy s, bool copy, bool ptr) {
  double height1 = 1.0;
  double height2 = 2.0;
  bool ok;

  // There is a huge amount of duplication here, but it's largely
  // unavoidable without missing out on doing what I want to do.
  if ( ptr ) {
    Plant p1(&s);
    p1.set_height(height1);

    if ( copy ) {
      Plant p2(p1);
      ok = p1 == p2;
    } else {
      Plant p2(&s);
      p2.set_height(height2);
      p2 = p1;
      ok = p1 == p2;
    }
  } else {
    Plant p1(s);
    p1.set_height(height1);

    if ( copy ) {
      Plant p2(p1);
      ok = p1 == p2;
    } else {
      Plant p2(s);
      p2.set_height(height2);
      p2 = p1;
      ok = p1 == p2;
    }
  }
  return ok;
}

}

}
