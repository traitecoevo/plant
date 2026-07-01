// -*-c++-*-
#ifndef PLANT_PLANT_INTERNALS_MINIMAL_H_
#define PLANT_PLANT_INTERNALS_MINIMAL_H_

#define HEIGHT_INDEX 0
#define MORTALITY_INDEX 1
#define FECUNDITY_INDEX 2

#include <memory> // std::shared_ptr
#include <odelia/ode_interface.hpp>
#include <vector>
// #include <plant/plant_internals.h>

// TODO(#483): extra_state bounds, upper and lower limits
namespace plant {

// Templated on the scalar type S so the per-plant state container can hold an
// AD active type for reverse-mode gradients (#472 scope B / #537, Milestone C).
// The `Internals` alias below pins S = double, so every existing use across the
// package keeps compiling and stays bit-identical; only AD paths instantiate
// basic_internals<ad_type>. NA_REAL initialisers are wrapped in S(...) so they
// also work when S is an AD type.
template <typename S>
class basic_internals {
public:
  basic_internals(size_t s_size=0, size_t a_size=0, size_t r_size=0)
      :
      state_size(s_size),
      aux_size(a_size),
      resource_size(r_size),
      states(s_size, S(0.0)),
      rates(s_size, S(NA_REAL)) ,
      auxs(a_size, S(0.0)),
      consumption_rates(r_size, S(NA_REAL))
    {}
  size_t state_size;
  size_t aux_size;
  size_t resource_size;


  // Perhaps make these private so the () overloads below have some use
  std::vector<S> states;
  std::vector<S> rates;
  std::vector<S> auxs;
  std::vector<S> consumption_rates;  // not quite as pithy

  S state(int i) const { return states[i]; }
  S rate(int i) const { return rates[i]; }
  S aux(int i) const { return auxs[i]; }
  S consumption_rate(int i) const { return consumption_rates[i]; }

  void set_state(int i, S v) { states[i] = v; }
  void set_rate(int i, S v) { rates[i] = v; }
  void set_aux(int i, S v) { auxs[i] = v; }
  void set_consumption_rate(int i, S v) { consumption_rates[i] = v; }

  void resize(size_t new_size, size_t new_aux_size) {
    state_size = new_size;
    aux_size = new_aux_size;
    states.resize(new_size, S(0.0));
    rates.resize(new_size, S(NA_REAL));
    auxs.resize(new_aux_size, S(0.0));
  }

  void resize_consumption_rates(size_t new_resource_size) {
    resource_size = new_resource_size;
    consumption_rates.resize(new_resource_size, S(NA_REAL));
  }
};

// Default state container used everywhere in the package (bit-identical to the
// pre-templating concrete class).
using Internals = basic_internals<double>;

} // namespace plant

#endif
