// -*-c++-*-
#ifndef PLANT_PLANT_CENSUS_H_
#define PLANT_PLANT_CENSUS_H_

#include <concepts>
#include <cstddef>
#include <string>
#include <vector>
#include <plant/internals.h>

namespace plant {

// One census metric: the per-individual quantity a census sums over the size
// distribution, and the name it is reported under. The quadrature is not here --
// it belongs to the reduction that owns the grid the density is carried in.
//
// The quantity is read from the strategy and one individual's own slots, which
// is the same pair every rate function reads: a strategy knows which slot holds
// what and reads it by index, where a caller outside would resolve a name.
//
// ⚠️ `of` MUST DECLARE ITS RETURN TYPE. An operator on an active scalar returns
// an expression template holding references to its operands, so a deduced return
// type hands back references to temporaries that are dead when the caller reads
// them -- and the sweep then walks whatever the reused stack holds. A plain
// function pointer is what makes that a compile error here rather than a
// conversion that succeeds.
template <class Strategy>
struct census_metric {
  using value_type = typename Strategy::value_type;
  const char* name;
  value_type (*of)(const Strategy&, const Internals<value_type>&);
};

// A strategy that says what a census of it reads. Asserted where a census is
// taken rather than here: a model nothing censuses owes this nothing, and a
// model that is censused says so at the call site with the member named instead
// of failing inside whatever loop reached for it.
template <class Strategy>
concept Censusable = requires {
  { *Strategy::census_metrics().begin() } ->
      std::convertible_to<census_metric<Strategy>>;
};

// The names of a strategy's metrics, in the order it reports them. Read off the
// metrics themselves, so a name and the quantity it belongs to cannot be listed
// in two places.
template <class Strategy>
  requires Censusable<Strategy>
std::vector<std::string> census_metric_names() {
  std::vector<std::string> ret;
  for (const census_metric<Strategy>& m : Strategy::census_metrics()) {
    ret.push_back(m.name);
  }
  return ret;
}

}  // namespace plant

#endif
