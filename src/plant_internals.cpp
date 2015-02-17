#include <tree2/plant_internals.h>
#include <R.h>

namespace tree2 {

Plant_internals::Plant_internals()
  : leaf_mass(NA_REAL),
    leaf_area(NA_REAL),
    height(NA_REAL),
    sapwood_mass(NA_REAL),
    bark_mass(NA_REAL),
    heartwood_area(0),
    heartwood_mass(0),
    root_mass(NA_REAL),
    live_mass(NA_REAL),
    assimilation(NA_REAL),
    respiration(NA_REAL),
    turnover(NA_REAL),
    net_production(NA_REAL),
    reproduction_fraction(NA_REAL),
    fecundity_rate(NA_REAL),
    leaf_area_growth_rate(NA_REAL),
    height_growth_rate(NA_REAL),
    heartwood_area_rate(NA_REAL),
    heartwood_mass_rate(NA_REAL),
    mortality_rate(NA_REAL),
    // But these should be zero
    mortality(0.0),
    fecundity(0.0) {
}

}
