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
template <class M>
struct census_metric {
  using value_type = typename M::value_type;
  const char* name;
  value_type (*of)(const M&, const Internals<value_type>&);
};

// A strategy that says what a census of it reads. Asserted where a census is
// taken rather than here: a model nothing censuses owes this nothing, and a
// model that is censused says so at the call site with the member named instead
// of failing inside whatever loop reached for it.
template <class M>
concept Censusable = requires {
  { *M::census_metrics().begin() } ->
      std::convertible_to<census_metric<M>>;
};

// A strategy whose fitness -- each species' net reproduction ratio and offspring
// production, read off the fecundity state -- is offered as a differentiable
// stand metric beside the census. Opt-in by a flag rather than inferred from the
// state, because a model has to say its fecundity channel is refereed: the
// reverse pass sweeps every state either way, and a row nobody checked would
// read as an answer.
template <class M>
concept DeclaresFitnessMetrics = requires {
  { M::differentiable_fitness } -> std::convertible_to<bool>;
} && bool(M::differentiable_fitness);

// The names of a strategy's metrics, in the order it reports them. Read off the
// metrics themselves, so a name and the quantity it belongs to cannot be listed
// in two places.
template <class M>
  requires Censusable<M>
std::vector<std::string> census_metric_names() {
  std::vector<std::string> ret;
  for (const census_metric<M>& m : M::census_metrics()) {
    ret.push_back(m.name);
  }
  return ret;
}

}  // namespace plant

#endif
