#include <tree2/plant.h>
#include <Rcpp.h>
#include <functional>

namespace tree2 {

Plant::Plant(strategy_type::ptr s)
  : strategy(s) {
  set_height(strategy->height_0);
}

// Individual size
// [eqn 1-8] Update size variables given an input leaf mass

// Height is really important -- everything else follows from it.
double Plant::height() const {
  return vars.height;
}
double Plant::height_dt() const {
  return vars.height_dt;
}
// NOTE: Only recomputes size variables if the height is actually
// different.  This is totally safe if nothing else sets either height
// or any size variable except this method.  This could help save
// quite a bit of calculation time and book-keeping down the track.
// If we're off because of a floating point difference, the worst that
// happens is that we recompute the variables again.
void Plant::set_height(double height_) {
  if (height_ < 0.0) {
    Rcpp::stop("height must be positive (given " +
               util::to_string(height_) + ")");
  }
  // TODO: in the original version of the model, all of plant size was driven
  // only  by height, so we only needed to check this here.  But now we have
  // two size variables (heartwood mass being the other) we need to be more
  // careful. For the new plant type, replace this height function with a new
  // function that takes both size variables as arguments and update all the
  // size variables at that point.

  if (!util::identical(height_, height())) {
    compute_vars_size(height_);
  }
}

double Plant::mortality() const {
  return vars.mortality;
}
double Plant::mortality_dt() const {
  return vars.mortality_dt;
}
void Plant::set_mortality(double x) {
  vars.mortality = x;
}

double Plant::fecundity() const {
  return vars.fecundity;
}
double Plant::fecundity_dt() const {
  return vars.fecundity_dt;
}
void Plant::set_fecundity(double x) {
  vars.fecundity = x;
}

double Plant::area_heartwood() const {
  return vars.area_heartwood;
}

double Plant::area_heartwood_dt() const {
  return vars.area_heartwood_dt;
}

void Plant::set_area_heartwood(double x) {
  // TODO: Consider recomputing the size variables here
  // See set_mass_heartwood
  vars.area_heartwood = x;
}

double Plant::mass_heartwood() const {
  return vars.mass_heartwood;
}

double Plant::mass_heartwood_dt() const {
  return vars.mass_heartwood_dt;
}

void Plant::set_mass_heartwood(double x) {
  // TODO: This needs to update size variables but does not yet,
  // and won't until we get the new version done.
  // See notes in set_height.
  vars.mass_heartwood = x;
}

// * Competitive environment
double Plant::area_leaf() const {
  return vars.area_leaf;
}

// [      ] Leaf area (not fraction) above height `z`
double Plant::area_leaf_above(double z) const {
  if (z < 0.0) {
    Rcpp::stop("Negative heights do not make sense");
  }
  return strategy->area_leaf_above(z, vars.height, vars.area_leaf);
}

// * Mass production
// [eqn 12-19,21] Update physiological variables given the current
// light environment (and given the current set of size variables).
void Plant::compute_vars_phys(const Environment& environment,
			      bool reuse_intervals) {
  // [eqn 12] Gross annual CO2 assimilation
  vars.assimilation = strategy->assimilation(environment, vars.height,
                                             vars.area_leaf, reuse_intervals);

  // [eqn 13] Total maintenance respiration
  vars.respiration = strategy->respiration(vars.mass_leaf,
                                           vars.mass_sapwood,
                                           vars.mass_bark,
                                           vars.mass_root);

  // [eqn 14] Total turnover
  vars.turnover = strategy->turnover(vars.mass_leaf, vars.mass_sapwood,
                                     vars.mass_bark, vars.mass_root);

  // [eqn 15] Net production
  vars.net_mass_production_dt = strategy->net_mass_production_dt(vars.assimilation,
                                                 vars.respiration,
                                                 vars.turnover);

  if (vars.net_mass_production_dt > 0) {
    // [eqn 16] - Fraction of whole plant growth that is leaf
    vars.fraction_allocation_reproduction =
      strategy->fraction_allocation_reproduction(vars.height);

    // [eqn 17] - Rate of offspring production
    //
    // NOTE: In EBT, was multiplied by Pi_0 (survival during
    // dispersal), but we do not here.
    // NOTE: This is also a hyperparametrisation and should move into
    // the initialisation function.
    vars.fecundity_dt =
      strategy->fecundity_dt(vars.net_mass_production_dt,
                              vars.fraction_allocation_reproduction);

    // [eqn 19] - Growth rate in leaf height
    // different to Falster 2010, which was growth rate in leaf mass
    vars.darea_leaf_dmass_live =
      strategy->darea_leaf_dmass_live(vars.area_leaf);
    vars.fraction_allocation_growth = strategy->fraction_allocation_growth(vars.height);

    vars.area_leaf_dt =
      vars.net_mass_production_dt * vars.fraction_allocation_growth *
      vars.darea_leaf_dmass_live;
   vars.height_dt =
      strategy->dheight_darea_leaf(vars.area_leaf) *
      vars.area_leaf_dt;
   vars.area_heartwood_dt =
     strategy->area_heartwood_dt(vars.area_leaf);
   vars.mass_heartwood_dt =
      strategy->mass_heartwood_dt(vars.mass_sapwood);
  } else {
    vars.fraction_allocation_reproduction = 0.0;
    vars.fraction_allocation_growth       = 0.0;
    vars.fecundity_dt        = 0.0;
    vars.area_leaf_dt = 0.0;
    vars.darea_leaf_dmass_live = 0.0;
    vars.height_dt    = 0.0;
    vars.area_heartwood_dt   = 0.0;
    vars.mass_heartwood_dt   = 0.0;
  }

  // [eqn 21] - Instantaneous mortality rate
  vars.mortality_dt =
      strategy->mortality_dt(vars.net_mass_production_dt / vars.area_leaf,
                             vars.mortality);
}

// Extra accounting.
// TODO: This will move into the "super size" plant.
void Plant::compute_vars_growth() {
  const strategy_type *s = strategy.get(); // for brevity.
  const double area_leaf = vars.area_leaf,
    area_leaf_dt = vars.area_leaf_dt;

  const double area_stem = vars.area_bark +
                           vars.area_sapwood +
                           vars.area_heartwood;

  // Changes with leaf area:
  vars.dheight_darea_leaf       = s->dheight_darea_leaf(area_leaf);
  vars.dmass_sapwood_darea_leaf = s->dmass_sapwood_darea_leaf(area_leaf);
  vars.dmass_bark_darea_leaf    = s->dmass_bark_darea_leaf(area_leaf);
  vars.dmass_root_darea_leaf    = s->dmass_root_darea_leaf(area_leaf);

  // Changes over time:
  vars.area_sapwood_dt    = s->area_sapwood_dt(area_leaf_dt);
  vars.area_bark_dt       = s->area_bark_dt(area_leaf_dt);
  vars.area_stem_dt       = s->area_stem_dt(area_leaf,
                                              area_leaf_dt);
  vars.diameter_stem_dt   = s->diameter_stem_dt(area_stem,
                                                  vars.area_stem_dt);
  vars.mass_root_dt       = s->mass_root_dt(area_leaf,
                                             area_leaf_dt);
  vars.mass_live_dt       = s->mass_live_dt(vars.fraction_allocation_reproduction,
                                              vars.net_mass_production_dt);
  vars.mass_total_dt      = s->mass_total_dt(vars.fraction_allocation_reproduction,
                                              vars.net_mass_production_dt,
                                              vars.mass_heartwood_dt);
  vars.mass_above_ground_dt =
    s->mass_above_ground_dt(area_leaf,
                             vars.fraction_allocation_reproduction,
                             vars.net_mass_production_dt,
                             vars.mass_heartwood_dt,
                             area_leaf_dt);

  // Odd one out:
  vars.ddiameter_stem_darea_stem = s->ddiameter_stem_darea_stem(area_stem);
}

// [eqn 20] Survival of seedlings during germination
double Plant::germination_probability(const Environment& environment) {
  return strategy->germination_probability(environment);
}

// ODE interface -- note that the don't care about time in the plant;
// only Patch and above does.
ode::const_iterator Plant::set_ode_state(ode::const_iterator it) {
  set_height(*it++);
  set_mortality(*it++);
  set_fecundity(*it++);
  set_area_heartwood(*it++);
  set_mass_heartwood(*it++);
  return it;
}
ode::iterator Plant::ode_state(ode::iterator it) const {
  *it++ = height();
  *it++ = mortality();
  *it++ = fecundity();
  *it++ = area_heartwood();
  *it++ = mass_heartwood();
  return it;
}
ode::iterator Plant::ode_rates(ode::iterator it) const {
  *it++ = height_dt();
  *it++ = mortality_dt();
  *it++ = fecundity_dt();
  *it++ = area_heartwood_dt();
  *it++ = mass_heartwood_dt();
  return it;
}

std::vector<std::string> Plant::ode_names() {
  return std::vector<std::string>({"height", "mortality", "fecundity",
	"area_heartwood", "mass_heartwood"});
}

// * R interface
Plant::strategy_type Plant::r_get_strategy() const {
  return *strategy.get();
}

Plant::Internals Plant::r_internals() const {
  return vars;
}

// This exists only so that I know that nothing will change the
// control parameters by only being able to access a const reference
// (it's shared with everything else that shares the strategy).  It
// also saves a little ugly looking referencing.
const Control& Plant::control() const {
  return strategy->control;
}

// * Private methods

// * Individual size
// [eqn 1-8] Update size variables to a new leaf mass.
void Plant::compute_vars_size(double height_) {
  vars.height = height_;
  vars.area_leaf = strategy->area_leaf(vars.height);
  vars.mass_leaf = strategy->mass_leaf(vars.area_leaf);
  vars.area_sapwood =  strategy->area_sapwood(vars.area_leaf);
  vars.mass_sapwood =  strategy->mass_sapwood(vars.area_sapwood, vars.height);
  vars.area_bark =  strategy->area_bark(vars.area_leaf);
  vars.mass_bark = strategy->mass_bark(vars.area_bark, vars.height);
  vars.mass_root = strategy->mass_root(vars.area_leaf);
  vars.mass_live = strategy->mass_live(vars.mass_leaf, vars.mass_sapwood,
                                       vars.mass_bark , vars.mass_root);
  vars.area_stem = strategy->area_stem(vars.area_bark, vars.area_sapwood,
                                         vars.area_heartwood);
  vars.mass_total = strategy->mass_total(vars.mass_leaf, vars.mass_bark, vars.mass_sapwood,
                                         vars.mass_heartwood, vars.mass_root);
  vars.mass_above_ground = strategy->mass_above_ground(vars.mass_leaf, vars.mass_bark,
                                                       vars.mass_sapwood, vars.mass_root);
  vars.diameter_stem = strategy->diameter_stem(vars.area_stem);
}

Plant make_plant(Plant::strategy_type s) {
  return Plant(make_strategy_ptr(s));
}

}
