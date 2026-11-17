// -*-c++-*-
#ifndef PLANT_PLANT_PARS_TABLE_H_
#define PLANT_PLANT_PARS_TABLE_H_

#include <array>
#include <cstddef>
#include <string>
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
// given a silent zero, so asking for one errors in these words. `read` is false
// for a parameter no equation reads: a trade-off may derive one and the chain
// rule loses nothing by its missing column, which it may not do for any other.
struct no_gradient {
  const char* name;
  const char* why;
  bool read = true;
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
// call it would be a string comparison against every excluded name for each
// field, and the gradient asks for this list once per recording.
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

// What a model's gradient hooks return, written once for every model: the
// addresses and names of P's columns, then the birth rate's where it is one
// (`birth_rate` null where it is not), and P's excluded parameters.
template <class P, class S>
std::vector<S*> column_ptrs(P& pars, S* birth_rate) {
  std::vector<S*> ret;
  ret.reserve(P::column_count + (birth_rate != nullptr ? 1 : 0));
  for (const auto& p : P::ad_parameter_columns) {
    ret.push_back(&(pars.*p.at));
  }
  if (birth_rate != nullptr) {
    ret.push_back(birth_rate);
  }
  return ret;
}

template <class P>
std::vector<std::string> column_names(bool birth_rate) {
  std::vector<std::string> ret;
  ret.reserve(P::column_count + (birth_rate ? 1 : 0));
  for (const auto& p : P::ad_parameter_columns) {
    ret.push_back(p.name);
  }
  if (birth_rate) {
    ret.push_back("birth_rate");
  }
  return ret;
}

template <class P>
std::vector<no_gradient> undifferentiable() {
  return std::vector<no_gradient>(P::undifferentiable.begin(),
                                  P::undifferentiable.end());
}

}  // namespace pars_table
}  // namespace plant

#endif
