// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

// Compatibility shim. The leaf gas-exchange and hydraulics model used to live
// here (inst/include/plant/leaf_model.h + src/leaf_model.cpp); it now ships as
// the standalone, header-only `phylloptim` package and plant consumes it via
// LinkingTo. See that package's README and PLAN.md.
//
// This header exists so that plant's own sources -- and in particular the ~17k
// lines of generated RcppR6/RcppExports glue, which name `plant::Leaf`
// throughout -- keep compiling unchanged. Aliasing costs nothing at runtime.
//
// Note that LinkingTo is NOT transitive in R: plant must name BH and odelia in
// its own DESCRIPTION even though it is `phylloptim` that includes them. It
// already does.

#include <phylloptim.hpp>

namespace plant {

// The leaf model itself.
using Leaf = ::phylloptim::Leaf;

// Constants that plant's own sources read out of the leaf model. Only one is
// actually used today (kg_per_mol_h2o, in TF24_Strategy::compute_rates, for the
// molar -> mass conversion on the soil water consumption rate), but the rest are
// pulled in so that the constants remain reachable as `plant::<name>` for
// anything that reaches for them later.
//
// Deliberately NOT re-exported: `phylloptim::gas_constant`, which was spelled `R` at
// plant namespace scope in the old header. A one-letter `R` in a public header
// is a collision hazard in a project where R is also the language and R_ prefixes
// its C API, and nothing outside the leaf model ever used it.
using ::phylloptim::C_to_K;
using ::phylloptim::gravity_head;
using ::phylloptim::H2O_CO2_stom_diff_ratio;
using ::phylloptim::kg_per_mol_h2o;
using ::phylloptim::kg_to_mol_h2o;
using ::phylloptim::kPa_to_Pa;
using ::phylloptim::umol_to_mol;

// Deliberately NOT re-exported: `umol_per_mol_to_Pa`. It was a namespace-scope
// constant 0.1013 = 1e-6 * 101300 Pa -- i.e. atmospheric pressure of 101.3 kPa in
// disguise, sitting next to a live, settable `atm_kpa_`. The leaf package derives
// it per-call as the member `Leaf::umol_per_mol_to_Pa_ = atm_kpa_ * kPa_to_Pa *
// umol_to_mol` (phylloptim item 10c), so there is no longer a namespace-scope
// constant to alias. Nothing in plant read it outside the leaf model.
//
// ⚠️ Because the conversion is now derived rather than fixed, `atm_kpa` reaches
// Gamma*, Kc, Ko, Km and the ci root-find's lower bound, which it did not before.
// TF24_Environment's driver therefore defaults to 101.3 (tf24_environment.h),
// the pressure the old constant encoded; it read 100.5 until this branch pinned
// it, which was worth ~0.8% relative on every TF24 run. See the NEWS entry.
// Setting the driver for altitude now moves those quantities too, as it should.

// Penman-Monteith leaf energy balance (#523).
using ::phylloptim::leaf_temp_max;
using ::phylloptim::leaf_temp_min;

} // namespace plant

#endif
