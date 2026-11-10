// -*-c++-*-
#ifndef PLANT_PLANT_PARS_TABLE_H_
#define PLANT_PLANT_PARS_TABLE_H_

#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

// What every model's parameter struct derives from its two lists. A struct P
// lists each of its members once, in `ad_parameter_fields`, and the members no
// gradient exists for, with the sentence that says why, in `undifferentiable`.
// Everything here is read off those two lists, so a member reaches all of it or
// none.
namespace plant {

// One registered parameter of P: where it lives and what its column is called.
// A member pointer rather than an address in one instance, which is what lets
// the table be constexpr.
template <class P, class S>
struct pars_field {
  S P::* at;
  const char* name;
};

// A parameter with no gradient, and why. It is excluded by name rather than
// given a silent zero, so asking for one errors in these words.
struct no_gradient {
  const char* name;
  const char* why;
};

namespace pars_table {

// A parameter carries a column unless `excluded` names it. Derived rather than
// declared: adding a member gives it a column, and taking one away is a
// deliberate entry with a sentence attached.
template <std::size_t N>
constexpr bool has_column(std::string_view name,
                          const std::array<no_gradient, N>& excluded) {
  for (const no_gradient& n : excluded) {
    if (name == n.name) {
      return false;
    }
  }
  return true;
}

// How many fields carry a gradient column, counted at compile time so a caller
// wanting the width does not build a vector of pointers to measure it.
template <class Fields, std::size_t N>
constexpr std::size_t column_count(const Fields& fields,
                                   const std::array<no_gradient, N>& excluded) {
  std::size_t n = 0;
  for (const auto& f : fields) {
    if (has_column(f.name, excluded)) ++n;
  }
  return n;
}

// The fields that carry a column, filtered once at compile time; filtered per
// call it was a string comparison against every excluded name for each field,
// and the gradient asks for this list once per recording.
template <std::size_t M, class Fields, std::size_t N>
constexpr std::array<typename Fields::value_type, M>
columns(const Fields& fields, const std::array<no_gradient, N>& excluded) {
  std::array<typename Fields::value_type, M> ret{};
  std::size_t n = 0;
  for (const auto& f : fields) {
    if (has_column(f.name, excluded)) {
      ret[n++] = f;
    }
  }
  return ret;
}

// A pointer to every member of `p`, in table order: the whole parameter set, so
// a rebind carries all of it across a scalar change.
template <class P, class S, std::size_t K>
std::vector<S*> field_ptrs(P& p, const std::array<pars_field<P, S>, K>& fields) {
  std::vector<S*> ret;
  ret.reserve(fields.size());
  for (const auto& f : fields) {
    ret.push_back(&(p.*f.at));
  }
  return ret;
}

}  // namespace pars_table
}  // namespace plant

#endif
