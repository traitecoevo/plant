#include <tree2/plant2.h>

namespace tree2 {

void Plant2::compute_vars_phys(const Environment& environment,
                               bool reuse_intervals) {
  // vars.height and vars.leaf_area are done on the way in here.  The
  // idea is to set height_rate, mortality_rate and fecundity_rate.
  const double height    = vars.height;
  const double leaf_area = vars.leaf_area;
  const double net_production =
    strategy->net_production(environment, height, leaf_area,
                               reuse_intervals);

  if (net_production > 0) {
    const double reproduction_fraction =
      strategy->reproduction_fraction(height);
    const double leaf_area_deployment_mass =
      strategy->leaf_area_deployment_mass(leaf_area);
    const double growth_fraction = strategy->growth_fraction(height);
    const double leaf_area_growth_rate =
      net_production * growth_fraction * leaf_area_deployment_mass;
    vars.height_rate =
      strategy->dheight_dleaf_area(leaf_area) * leaf_area_growth_rate;
    vars.fecundity_rate =
      strategy->dfecundity_dt(net_production,
                              reproduction_fraction);

  } else {
    vars.height_rate    = 0.0;
    vars.fecundity_rate = 0.0;
  }

  // [eqn 21] - Instantaneous mortality rate
  // TODO: Should this not be in Strategy @dfalster?  Move entire
  // if/else clause into Strategy, and rewrite to take net_production,
  // leaf_area and mortality as arguments.
  //
  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (R_FINITE(vars.mortality)) {
    vars.mortality_rate =
      strategy->mortality_dt(net_production / leaf_area);
  } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    vars.mortality_rate = 0.0;
  }
}

ode::const_iterator Plant2::set_ode_state(ode::const_iterator it) {
  set_height(*it++);
  set_mortality(*it++);
  set_fecundity(*it++);
  return it;
}
ode::iterator Plant2::ode_state(ode::iterator it) const {
  *it++ = height();
  *it++ = mortality();
  *it++ = fecundity();
  return it;
}
ode::iterator Plant2::ode_rates(ode::iterator it) const {
  *it++ = height_rate();
  *it++ = mortality_rate();
  *it++ = fecundity_rate();
  return it;
}

Plant2 make_plant2(Plant2::strategy_type s) {
  return Plant2(make_strategy_ptr(s));
}

}
