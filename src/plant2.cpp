#include <tree2/plant2.h>

namespace tree2 {

void Plant2::compute_vars_phys(const Environment& environment,
                               bool reuse_intervals) {
  // vars.height and vars.area_leaf are done on the way in here.  The
  // idea is to set height_dt, mortality_dt and reproduction_dt.
  const double height    = vars.height;
  const double area_leaf = vars.area_leaf;
  const double net_mass_production_dt =
    strategy->net_mass_production_dt(environment, height, area_leaf,
                               reuse_intervals);

  if (net_mass_production_dt > 0) {
    const double fraction_allocation_reproduction =
      strategy->fraction_allocation_reproduction(height);
    const double darea_leaf_dmass_live =
      strategy->darea_leaf_dmass_live(area_leaf);
    const double fraction_allocation_growth = strategy->fraction_allocation_growth(height);
    const double area_leaf_dt =
      net_mass_production_dt * fraction_allocation_growth * darea_leaf_dmass_live;
    vars.height_dt =
      strategy->dheight_darea_leaf(area_leaf) * area_leaf_dt;
    vars.reproduction_dt =
      strategy->reproduction_dt(net_mass_production_dt,
                              fraction_allocation_reproduction);

  } else {
    vars.height_dt    = 0.0;
    vars.reproduction_dt = 0.0;
  }

  // [eqn 21] - Instantaneous mortality rate

  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (R_FINITE(vars.mortality)) {
    vars.mortality_dt =
      strategy->mortality_dt(net_mass_production_dt / area_leaf);
  } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    vars.mortality_dt = 0.0;
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
  *it++ = height_dt();
  *it++ = mortality_dt();
  *it++ = reproduction_dt();
  return it;
}

Plant2 make_plant2(Plant2::strategy_type s) {
  return Plant2(make_strategy_ptr(s));
}

}
