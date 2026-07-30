// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

// Compatibility shim. The leaf gas-exchange and hydraulics model used to live
// here (inst/include/plant/leaf_model.h + src/leaf_model.cpp); it now ships as
// the standalone, header-only `leaf` package and plant consumes it via
// LinkingTo. See that package's README and PLAN.md.
//
// This header exists so that plant's own sources -- and in particular the ~17k
// lines of generated RcppR6/RcppExports glue, which name `plant::Leaf`
// throughout -- keep compiling unchanged. Aliasing costs nothing at runtime.
//
// Note that LinkingTo is NOT transitive in R: plant must name BH and odelia in
// its own DESCRIPTION even though it is `leaf` that includes them. It already
// does.

#include <leaf.hpp>

namespace plant {

// The leaf model itself.
using Leaf = ::leaf::Leaf;

// Constants that plant's own sources read out of the leaf model. Only one is
// actually used today (kg_per_mol_h2o, in TF24_Strategy::compute_rates, for the
// molar -> mass conversion on the soil water consumption rate), but the rest are
// pulled in so that the constants remain reachable as `plant::<name>` for
// anything that reaches for them later.
//
// Deliberately NOT re-exported: `leaf::gas_constant`, which was spelled `R` at
// plant namespace scope in the old header. A one-letter `R` in a public header
// is a collision hazard in a project where R is also the language and R_ prefixes
// its C API, and nothing outside the leaf model ever used it.
using ::leaf::C_to_K;
using ::leaf::gravity_head;
using ::leaf::H2O_CO2_stom_diff_ratio;
using ::leaf::kg_per_mol_h2o;
using ::leaf::kg_to_mol_h2o;
using ::leaf::kPa_to_Pa;
using ::leaf::umol_per_mol_to_Pa;
using ::leaf::umol_to_mol;

// Penman-Monteith leaf energy balance (#523).
using ::leaf::leaf_temp_max;
using ::leaf::leaf_temp_min;

} // namespace plant

#endif
