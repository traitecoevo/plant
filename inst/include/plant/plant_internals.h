// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_INTERNALS_H_
#define PLANT_PLANT_PLANT_INTERNALS_H_

namespace plant {

// These are common to all minimal plants, for now at least.
//
// Moving to a more general "size" based model would be easy enough
// but we'd need to also store height because Patch & Environment
// between them use height to work out how far up to compute the
// canopy openness for.  So like leaf_area being carried around we'd
// need to carry height as well.
struct Plant_internals {
  Plant_internals()
    :
    height(NA_REAL),
    height_dt(NA_REAL),
    mortality(0.0),
    mortality_dt(NA_REAL),
    fecundity(0.0),
    fecundity_dt(NA_REAL) {
  }
  double height;
  double area_leaf;
  double height_dt;
  double mortality;
  double mortality_dt;
  double fecundity;
  double fecundity_dt;
};

}

#endif
