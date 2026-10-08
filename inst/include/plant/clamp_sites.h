// -*-c++-*-
#ifndef PLANT_PLANT_CLAMP_SITES_H_
#define PLANT_PLANT_CLAMP_SITES_H_

#include <phylloptim/clamp_sites.hpp>
#include <cstddef>
#include <memory>
#include <vector>

namespace plant {

// Where the model replaces a value by a bound. Each one holds the quantity it
// masks constant, so whatever the masked value fed stops responding -- and a row
// severed that way and a row that is honestly zero are the same number. Counting
// is what separates them, which is why a site is added to this list rather than
// given a counter of its own.
//
// The sites live in two objects and the list is one, so a tally reads the same
// however the site is reached.
enum clamp_site {
  // The light floor, at two depths of the same constant. The crown one floors
  // each quadrature point of the mean-light integrand and the radiation one
  // floors the mean those points make -- so on the shipped shading model the
  // crown site binds FIRST and the radiation site sees an already-floored mean.
  CLAMP_LIGHT_FLOOR = 0,
  CLAMP_LIGHT_FLOOR_CROWN,
  // The soil retention curve, which diverges as moisture goes to zero.
  CLAMP_SOIL_MOISTURE_FLOOR,
  CLAMP_SOIL_POTENTIAL_CEILING,
  // A layer held at zero rate because it is at the residual and still drying.
  CLAMP_SOIL_POSITIVITY,
  // Hydraulic conductivity, floored at zero and ceilinged at saturation.
  CLAMP_SOIL_CONDUCTIVITY,
  // The two water inputs, each a positive part over a spline that can undershoot.
  CLAMP_RAINFALL,
  CLAMP_INFILTRATION,
  // Rooting depth capped: where it binds, the root profile stops responding to
  // height, which is a state row rather than a numeric guard.
  CLAMP_ROOTING_DEPTH,
  // The leaf model's own sites, appended in phylloptim's enum order so that the two
  // lists fold by offset rather than by a name lookup. They are counted differently
  // from everything above: the leaf solves in double on BOTH paths, so which path a
  // clamp fired on is a question of WHEN rather than of the scalar, and the answer
  // comes from a delta taken across record_leaf_outputs.
  //
  // The block's width is the leaf model's own count rather than a number repeated
  // here: a site added there would otherwise land past CLAMP_SITE_COUNT, and every
  // reader that stops at the count drops it without saying so.
  CLAMP_LEAF_FIRST,
  CLAMP_SITE_COUNT = CLAMP_LEAF_FIRST + phylloptim::CLAMP_SITE_COUNT
};

inline const char* clamp_site_name(int site) {
  switch (site) {
  case CLAMP_LIGHT_FLOOR:            return "light_floor";
  case CLAMP_LIGHT_FLOOR_CROWN:      return "light_floor_crown";
  case CLAMP_SOIL_MOISTURE_FLOOR:    return "soil_moisture_floor";
  case CLAMP_SOIL_POTENTIAL_CEILING: return "soil_potential_ceiling";
  case CLAMP_SOIL_POSITIVITY:        return "soil_positivity";
  case CLAMP_SOIL_CONDUCTIVITY:      return "soil_conductivity";
  case CLAMP_RAINFALL:               return "rainfall";
  case CLAMP_INFILTRATION:           return "infiltration";
  case CLAMP_ROOTING_DEPTH:          return "rooting_depth";
  }
  // The leaf block's names are the leaf model's, taken by offset. Spelling them
  // again here would let the two lists disagree while both compiled.
  if (site >= CLAMP_LEAF_FIRST && site < CLAMP_SITE_COUNT) {
    return phylloptim::clamp_site_name(site - CLAMP_LEAF_FIRST);
  }
  return "unknown";
}

// One object's tally, kept as two because the two paths answer different
// questions: the forward one says the guard is reachable, and only the
// differentiated one says a gradient carries a severance.
//
// The differentiated half is held behind a pointer so a rebind can SHARE it. A
// block copies its strategy and environment per unit and discards them, so a
// plain member would take every count on that path down with the copy -- which
// is what made the site look unreachable while it was firing a million times.
struct clamp_counter {
  std::vector<std::size_t> forward =
    std::vector<std::size_t>(CLAMP_SITE_COUNT, 0);
  std::shared_ptr<std::vector<std::size_t>> differentiated =
    std::make_shared<std::vector<std::size_t>>(CLAMP_SITE_COUNT, 0);

  void clear() {
    forward.assign(CLAMP_SITE_COUNT, 0);
    differentiated->assign(CLAMP_SITE_COUNT, 0);
  }
};

}  // namespace plant

#endif
