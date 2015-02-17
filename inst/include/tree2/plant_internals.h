// -*-c++-*-
#ifndef TREE2_PLANT_INTERNALS_H_
#define TREE2_PLANT_INTERNALS_H_

namespace tree2 {

struct Plant_internals {
public:
  Plant_internals();
  // * Individual size
  // Mass of leaves.  This is the core independent variable
  double leaf_mass;      // [eqn 1]
  // Other size variables that follow directly from `leaf_mass`:
  double leaf_area;      // [eqn 2]
  double height;         // [eqn 3]
  double sapwood_area;   // [eqn 4]
  double sapwood_mass;   // [eqn 4]
  double bark_area;      // [eqn 5]
  double bark_mass;      // [eqn 5]
  double heartwood_area;
  double heartwood_mass; // [eqn 6]
  double basal_area;   //TODO: stem area?
  double root_mass;      // [eqn 7] (fine roots)
  double live_mass;      // [eqn 8]
  double total_mass;
  double above_ground_mass;
  double diameter;      //TODO: stem_diameter?

  // * Mass production
  double assimilation;   // [eqn 12] Gross annual CO2 assimilation
  double respiration;    // [eqn 13] Total maintenance respiration
  double turnover;       // [eqn 14] Total turnover
  double net_production; // [eqn 15] Net production
  double reproduction_fraction; // [eqn 16]
  double growth_fraction;
  double fecundity_rate; // [eqn 17] Rate of offspring production
  double leaf_area_growth_rate; // [eqn 19] Growth rate in leaf mass
  double leaf_area_deployment_mass;
  double height_growth_rate;    // [doc/details.md]
  double heartwood_area_rate;
  double heartwood_mass_rate;
  // * Mortality
  double mortality_rate; // [eqn 21]
  // * Variables
  double mortality;
  double fecundity;

  // * Growth partitioning
  double dheight_dleaf_area;
  double dsapwood_mass_dleaf_area;
  double dbark_mass_dleaf_area;
  double droot_mass_dleaf_area;
  double dsapwood_area_dt;
  double dbark_area_dt;
  double dbasal_area_dt;
  double dbasal_diam_dt;
  double droot_mass_dt;
  // TODO: reproduction_fraction -> reproduction_mass_fraction?
  // TODO: net production -> net mass production?
  double dlive_mass_dt;
  double dtotal_mass_dt;
  double dabove_ground_mass_dt;
  double dbasal_diam_dbasal_area;
};

}

#endif
