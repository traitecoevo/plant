// -*-c++-*-
#ifndef PLANT_PLANT_WITH_SLOPE_H_
#define PLANT_PLANT_WITH_SLOPE_H_

#include <odelia/value_with_slope.hpp>

namespace plant {

// ⚠️ THE TYPE IS ODELIA'S, AND IT HAS TO BE. Its `for_each_active` exists for
// odelia's own visit_active, which passes over any shape it does not open
// WITHOUT refusing it -- so a pair that does not say what it holds loses both
// members from the tape audit in silence, and active_system::release is what
// reports the miss. A model defining the pair itself is a model carrying that
// library's obligation, and a second definition is a second thing to forget.
template <typename T>
using with_slope = odelia::value_with_slope<T>;

// WHAT THE SLOPE IS WITH RESPECT TO is this package's, which is why odelia's
// header does not name it: phylloptim's operating point carries d/d(collar
// potential) in the same type.
//
// Here it is the vertical slope, which is what every frame of the competition
// path carries: an individual's contribution, a species' sum over its nodes, a
// patch's sum over its species, and a light-environment knot.
//
// `slope` is the SIGNED d(value)/d(height), so it is negative wherever the
// profile falls with height. The canopy shape forms the pair the other way round
// -- q is exactly -dQ/dz -- and each strategy's compute_competition_and_slope
// negates it on the way out, which is the one place the convention turns over.
// Above there it is signed, which is what build_extinction_field's chain rule
// assumes.
//
// The shape's own pair stays a pair of Q and q: those are the model's names for
// two quantities it takes from one u^eta, not a value and a slope.

}

#endif
