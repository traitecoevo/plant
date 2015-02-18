#include <tree2/plant2.h>

namespace tree2 {

void Plant2::compute_vars_phys(const Environment& environment,
                               bool reuse_intervals) {
  // vars.height and vars.area_leaf are done on the way in here.  The
  // idea is to set height_dt, mortality_dt and fecundity_dt.
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
    vars.fecundity_dt =
      strategy->fecundity_dt(net_mass_production_dt,
                              fraction_allocation_reproduction);

  } else {
    vars.height_dt    = 0.0;
    vars.fecundity_dt = 0.0;
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.mortality_dt =
      strategy->mortality_dt(net_mass_production_dt / area_leaf,
                             vars.mortality);
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
  *it++ = fecundity_dt();
  return it;
}

Plant2 make_plant2(Plant2::strategy_type s) {
  return Plant2(make_strategy_ptr(s));
}

}
