#include <tree/ffw16_plant_plus_internals.h>
#include <R.h>

namespace tree {

// TODO: Some of these start at zero, others at NA_REAL.  Work out
// when they get set and set them.
//
// TODO: Some are d<foo>_dt, others are foo_rate.  I'm OK with time
// derivatives getting separate treatment, but it needs to be
// consistent.
FFW16_PlantPlus_internals::FFW16_PlantPlus_internals()
  :
  // * Size
  mass_leaf(NA_REAL),
  area_leaf(NA_REAL),
  height(NA_REAL),
  area_sapwood(NA_REAL),
  mass_sapwood(NA_REAL),
  area_bark(NA_REAL),
  mass_bark(NA_REAL),
  area_heartwood(0.0), // NOTE: Starts at zero
  mass_heartwood(0.0), // NOTE: Starts at zero
  area_stem(NA_REAL),
  mass_root(NA_REAL),
  mass_live(NA_REAL),
  mass_total(NA_REAL),
  mass_above_ground(NA_REAL),
  diameter_stem(NA_REAL),
  // * Physiological
  assimilation(NA_REAL),
  respiration(NA_REAL),
  turnover(NA_REAL),
  net_mass_production_dt(NA_REAL),
  fraction_allocation_reproduction(NA_REAL),
  fraction_allocation_growth(NA_REAL),
  fecundity_dt(NA_REAL),
  area_leaf_dt(NA_REAL),
  darea_leaf_dmass_live(NA_REAL),
  height_dt(NA_REAL),
  area_heartwood_dt(NA_REAL),
  mass_heartwood_dt(NA_REAL),
  mortality_dt(NA_REAL),
  // But these should be zero
  mortality(0.0),
  fecundity(0.0),
  // * Growth
  dheight_darea_leaf(NA_REAL),
  dmass_sapwood_darea_leaf(NA_REAL),
  dmass_bark_darea_leaf(NA_REAL),
  dmass_root_darea_leaf(NA_REAL),
  area_sapwood_dt(NA_REAL),
  area_bark_dt(NA_REAL),
  area_stem_dt(NA_REAL),
  diameter_stem_dt(NA_REAL),
  mass_root_dt(NA_REAL),
  mass_live_dt(NA_REAL),
  mass_total_dt(NA_REAL),
  mass_above_ground_dt(NA_REAL),
  ddiameter_stem_darea_stem(NA_REAL) {
}

}
