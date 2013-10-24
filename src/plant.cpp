#include "plant.h"

#include "find_root.h"

namespace model {

PlantBase::~PlantBase() {
}

Plant::Plant(Strategy s)
  : strategy(s) {
  set_height(strategy->height_0);
}

Plant::Plant(Strategy *s)
  : strategy(s) {
  set_height(strategy->height_0);
}

Plant::internals::internals() 
  : mass_leaf(NA_REAL),
    leaf_area(NA_REAL),
    height(NA_REAL),
    mass_sapwood(NA_REAL),
    mass_bark(NA_REAL),
    mass_heartwood(NA_REAL),
    mass_root(NA_REAL),
    mass_total(NA_REAL),
    // NOTE: Zero to allow '==' to work on new plants
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
bool Plant::operator==(const Plant &rhs) const {
  if (strategy.standalone() && strategy.get() == rhs.strategy.get())
    ::Rf_warning("This is going to end badly.");
  return strategy == rhs.strategy && vars == rhs.vars;
}

// NOTE: Using identical to compare between floats here -- we *want*
// identical comparisons for the way this method is used (are objects
// absolutely identical?).
bool Plant::internals::operator==(const Plant::internals &rhs) const {
  const bool ret = true
    // Size members...
    && util::identical(mass_leaf,         rhs.mass_leaf)
    && util::identical(leaf_area,         rhs.leaf_area)
    && util::identical(height,            rhs.height)
    && util::identical(mass_sapwood,      rhs.mass_sapwood)
    && util::identical(mass_bark,         rhs.mass_bark)
    && util::identical(mass_heartwood,    rhs.mass_heartwood)
    && util::identical(mass_root,         rhs.mass_root)
    && util::identical(mass_total,        rhs.mass_total)
    // ...physiological members...
    && util::identical(assimilation,      rhs.assimilation)
    && util::identical(respiration,       rhs.respiration)
    && util::identical(turnover,          rhs.turnover)
    && util::identical(net_production,    rhs.net_production)
    && util::identical(reproduction_fraction,    rhs.reproduction_fraction)
    && util::identical(fecundity_rate,    rhs.fecundity_rate)
    && util::identical(leaf_fraction,     rhs.leaf_fraction)
    && util::identical(mass_leaf_growth_rate,    rhs.mass_leaf_growth_rate)
    && util::identical(height_growth_rate,       rhs.height_growth_rate)
    && util::identical(mortality_rate,    rhs.mortality_rate)
    // ...and "other" members...
    && util::identical(fecundity,         rhs.fecundity)
    && util::identical(mortality,         rhs.mortality)
    ;
  return ret;
}

double Plant::height() const {
  return vars.height;
}

// [eqn 1-8] Update size variables given an input leaf mass
//
// NOTE: Only recomputes size variables if the height is actually
// different.  This is totally safe if nothing else sets either height
// or any size variable except this method.  This could help save
// quite a bit of calculation time and book-keeping down the track.
// If we're off because of a floating point difference, the worst that
// happens is that we recompute the variables again.
void Plant::set_height(double height_) {
  if (height_ < 0.0)
    ::Rf_error("height must be positive (given %2.5f)", height_);
  if (!util::identical(height_, height()))
    compute_vars_size(height_);
}

double Plant::height_rate() const {
  return vars.height_growth_rate;
}

// * Competitive environment
// [      ] Leaf area (not fraction) above height `z`
double Plant::leaf_area_above(double z) const {
  if (z < 0.0)
    ::Rf_error("Negative heights do not make sense");
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

  if (vars.net_production > 0) {
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

    // [eqn 19] - Growth rate in leaf mass
    vars.mass_leaf_growth_rate = vars.net_production *
      (1 - vars.reproduction_fraction) * vars.leaf_fraction;

    // [      ] - see doc/details.md
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
  if (vars.fecundity > 1)
    vars.fecundity = modf(vars.fecundity, &born);
  return static_cast<int>(born);
}

bool Plant::died() {
  const int did_die = unif_rand() < mortality_probability();
  set_mortality(0.0);
  return did_die;
}

// [eqn 20] Survival of seedlings during germination
//
// NOTE: This does not check/enforce that height is set to the seed
// height (so this is actually the germination probability of a plant
// that happens to be the current size).  This might be something to
// change.
double Plant::germination_probability(const Environment& environment) {
  compute_vars_phys(environment);
  if (vars.net_production > 0) {
    const double tmp = vars.leaf_area * strategy->c_s0 / vars.net_production;
    return 1 / (tmp * tmp + 1.0);
  } else {
    return 0.0;
  }
}

// * ODE interface
size_t Plant::ode_size() const { 
  return ode_dimension; 
}

ode::iterator_const Plant::set_ode_values(double /* unused: time */,
					  ode::iterator_const it) {
  set_height(*it++);
  set_mortality(*it++);
  set_fecundity(*it++);
  return it;
}

ode::iterator Plant::ode_values(ode::iterator it) const {
  *it++ = height();
  *it++ = mortality();
  *it++ = fecundity();
  return it;
}

ode::iterator Plant::ode_rates(ode::iterator it) const {
  *it++ = height_rate();
  *it++ = mortality_rate();
  *it++ = fecundity_rate();
  return it;
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

// * Access the Control parameter (protected)
//
// This exists only so that I know that nothing will change the
// control parameters by only being able to access a const reference
// (it's shared with everything else that shares the strategy).  It
// also saves a little ugly looking referencing.
const Control& Plant::control() const {
  return strategy->control;
}

// This is all fairly unfortunate.  Perhaps it will change later?
integration::intervals_type Plant::get_last_integration_intervals() const {
  return strategy->integrator.get_last_intervals();
}

void Plant::set_integration_intervals(integration::intervals_type x) {
  integration_intervals = x;
}

// * Private methods

// * Individual size
// [eqn 1-8] Update size variables to a new leaf mass.
void Plant::compute_vars_size(double height_) {
  // First 3 differ from paper; working height->mass, not mass->height.
  // [eqn 3] height
  vars.height = height_;
  // [eqn 2] leaf_area (inverse of [eqn 3])
  vars.leaf_area = pow(vars.height / strategy->a1, 1 / strategy->B1);
  // [eqn 1] mass_leaf (inverse of [eqn 2])
  vars.mass_leaf = vars.leaf_area * strategy->lma;

  // These are identical to paper.
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
  if (z > vars.height)
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
// NOTE: In contrast with Daniel's implementation (but following
// Falster 2012), we do not normalise by Y*c_bio here.
double Plant::compute_assimilation(const Environment& environment) {
  FunctorBind1<Plant, const Environment&,
	       &Plant::compute_assimilation_x> fun(this, environment);
  const bool over_distribution =
    control().plant_assimilation_over_distribution;
  const double x_min = 0, x_max = over_distribution ? 1 : vars.height;
  double A = 0.0;

  // This whole section exists only because of working around
  // CohortTop, where we want to recompute the photosynthetic integral
  // using the same grid that we refined onto before.  This would
  // probably be heaps more efficient if we worked with references or
  // iterators and avoided the copying that is going on here.
  // Conceptually, this probably belongs in its own function.
  if (control().plant_assimilation_reuse_intervals &&
      integration_intervals.size() == 2) {
    if (over_distribution) {
      A = strategy->integrator.integrate_with_intervals(&fun,
							integration_intervals);
    } else {
      // In theory we might not have to to scale if we're doing
      // backward differencing.
      A = strategy->integrator.integrate_with_intervals(&fun,
							integration::internal::rescale_intervals(integration_intervals, x_min, x_max));
    }
  } else {
    A = strategy->integrator.integrate(&fun, x_min, x_max);
  }
  return vars.leaf_area * A;
}

// This is used in the calculation of assimilation by
// `compute_assimilation` above; it is the term within the integral in
// [eqn 12]; i.e., A_lf(A_0v, E(z,a)) * q(z,h(m_l))
// where `z` is height.
double Plant::compute_assimilation_x(double x, 
				     const Environment& environment) const {
  if (control().plant_assimilation_over_distribution)
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

// [eqn 18] Fraction of mass growth that is leaves (see doc/details.md
// for derivation).
double Plant::compute_leaf_fraction() const {
  const Strategy *s = strategy.get(); // for brevity.
  return 1.0/(// d m_l / d m_l
	      1.0 +
	      // d m_s / d m_l + d m_b / d m_l:
	      (1.0 + s->b) *
	      s->rho * s->eta_c * s->a1 / (s->theta * s->lma) *
	      (s->B1 + 1.0) * pow(vars.leaf_area, s->B1) +
	      // d m_h / d m_l:
	      s->rho * s->eta_c * s->a2 / vars.mass_leaf *
	      s->B2 * pow(vars.leaf_area, s->B2) +
	      // d m_r / d m_l:
	      s->a3 / s->lma);
}

// NOTE: static method
double Plant::height_seed(Strategy *s) {
  Plant p(s);
  const double mass_seed = s->s;
  const double
    h0 = p.height_given_mass_leaf(DBL_MIN),
    h1 = p.height_given_mass_leaf(mass_seed);

  // Functor computing total mass of a seed with height h, minus the
  // target seed mass.
  util::FunctorRoot<Plant, &Plant::mass_total_given_height>
    fun(&p, mass_seed);

  const double tol = p.control().plant_seed_tol;
  const int max_iterations = p.control().plant_seed_iterations;
  util::RootFinder root(tol, tol, max_iterations);

  return root.root(&fun, h0, h1);
}

// NOTE: This is used only by height_seed.
double Plant::mass_total_given_height(double h) {
  set_height(h);
  return vars.mass_total;
}

// NOTE: static method.
void Plant::prepare_strategy(Strategy *s) {
  s->eta_c = 1 - 2/(1 + s->eta) + 1/(1 + 2*s->eta);
  s->k_l = s->a4 * pow(s->lma, -s->B4);
  s->height_0 = height_seed(s);
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
			       _["height_growth_rate"]=
			       vars.height_growth_rate,
			       _["mortality_rate"]=vars.mortality_rate);
}

bool Plant::r_died() {
  Rcpp::RNGScope scope;
  return died();
}

// This is useful for finding the seed height.
double Plant::height_given_mass_leaf(double mass_leaf_) const {
  return strategy->a1 * pow(mass_leaf_ / strategy->lma, strategy->B1);
}

namespace test {

bool test_plant(Strategy s, bool copy, bool ptr) {
  double height1 = 1.0;
  double height2 = 2.0;
  bool ok;

  // There is a huge amount of duplication here, but it's largely
  // unavoidable without missing out on doing what I want to do.
  if (ptr) {
    Plant p1(&s);
    p1.set_height(height1);

    if (copy) {
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

    if (copy) {
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
