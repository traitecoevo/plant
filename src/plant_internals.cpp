#include <tree2/plant_internals.h>
#include <R.h>

namespace tree2 {

// TODO: Some of these start at zero, others at NA_REAL.  Work out
// when they get set and set them.
//
// TODO: Some are d<foo>_dt, others are foo_rate.  I'm OK with time
// derivatives getting separate treatment, but it needs to be
// consistent.
Plant_internals::Plant_internals()
  :
  // * Size
  leaf_mass(NA_REAL),
  leaf_area(NA_REAL),
  height(NA_REAL),
  sapwood_area(NA_REAL),
  sapwood_mass(NA_REAL),
  bark_area(NA_REAL),
  bark_mass(NA_REAL),
  heartwood_area(0.0), // NOTE: Starts at zero
  heartwood_mass(0.0), // NOTE: Starts at zero
  basal_area(NA_REAL),
  root_mass(NA_REAL),
  live_mass(NA_REAL),
  total_mass(NA_REAL),
  above_ground_mass(NA_REAL),
  diameter(NA_REAL),
  // * Physiological
  assimilation(NA_REAL),
  respiration(NA_REAL),
  turnover(NA_REAL),
  net_production(NA_REAL),
  reproduction_fraction(NA_REAL),
  growth_fraction(NA_REAL),
  fecundity_rate(NA_REAL),
  leaf_area_growth_rate(NA_REAL),
  leaf_area_deployment_mass(NA_REAL),
  height_growth_rate(NA_REAL),
  heartwood_area_rate(NA_REAL),
  heartwood_mass_rate(NA_REAL),
  mortality_rate(NA_REAL),
  // But these should be zero
  mortality(0.0),
  fecundity(0.0),
  // * Growth
  dheight_dleaf_area(NA_REAL),
  dsapwood_mass_dleaf_area(NA_REAL),
  dbark_mass_dleaf_area(NA_REAL),
  droot_mass_dleaf_area(NA_REAL),
  dsapwood_area_dt(NA_REAL),
  dbark_area_dt(NA_REAL),
  dbasal_area_dt(NA_REAL),
  dbasal_diam_dt(NA_REAL),
  droot_mass_dt(NA_REAL),
  dlive_mass_dt(NA_REAL),
  dtotal_mass_dt(NA_REAL),
  dabove_ground_mass_dt(NA_REAL),
  dbasal_diam_dbasal_area(NA_REAL) {
}

}
