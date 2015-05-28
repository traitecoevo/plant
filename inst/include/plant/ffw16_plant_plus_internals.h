// -*-c++-*-
#ifndef PLANT_PLANT_FFW16_PLANT_PLUS_INTERNALS_H_
#define PLANT_PLANT_FFW16_PLANT_PLUS_INTERNALS_H_

namespace plant {

struct FFW16_PlantPlus_internals {
public:
  FFW16_PlantPlus_internals();
  // * Individual size
  // Mass of leaves.  This is the core independent variable
  double mass_leaf;      // [eqn 1]
  // Other size variables that follow directly from `mass_leaf`:
  double area_leaf;      // [eqn 2]
  double height;         // [eqn 3]
  double area_sapwood;   // [eqn 4]
  double mass_sapwood;   // [eqn 4]
  double area_bark;      // [eqn 5]
  double mass_bark;      // [eqn 5]
  double area_heartwood;
  double mass_heartwood; // [eqn 6]
  double area_stem;
  double mass_root;      // [eqn 7] (fine roots)
  double mass_live;      // [eqn 8]
  double mass_total;
  double mass_above_ground;
  double diameter_stem;

  // * Mass production
  double assimilation;   // [eqn 12] Gross annual CO2 assimilation
  double respiration;    // [eqn 13] Total maintenance respiration
  double turnover;       // [eqn 14] Total turnover
  double net_mass_production_dt; // [eqn 15] Net production
  double fraction_allocation_reproduction; // [eqn 16]
  double fraction_allocation_growth;
  double fecundity_dt; // [eqn 17] Rate of offspring production
  double area_leaf_dt; // [eqn 19] Growth rate in leaf mass
  double darea_leaf_dmass_live;
  double height_dt;    // [doc/details.md]
  double area_heartwood_dt;
  double mass_heartwood_dt;
  // * Mortality
  double mortality_dt; // [eqn 21]
  // * Variables
  double mortality;
  double fecundity;

  // * Growth partitioning
  double dheight_darea_leaf;
  double dmass_sapwood_darea_leaf;
  double dmass_bark_darea_leaf;
  double dmass_root_darea_leaf;
  double area_sapwood_dt;
  double area_bark_dt;
  double area_stem_dt;
  double diameter_stem_dt;
  double mass_root_dt;
  double mass_live_dt;
  double mass_total_dt;
  double mass_above_ground_dt;
  double ddiameter_stem_darea_stem;
};

}

#endif
