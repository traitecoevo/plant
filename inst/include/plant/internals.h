// -*-c++-*-
#ifndef PLANT_PLANT_INTERNALS_MINIMAL_H_
#define PLANT_PLANT_INTERNALS_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <odelia/ode_interface.hpp>
#include <vector>
// #include <plant/plant_internals.h>

// TODO(#483): extra_state bounds, upper and lower limits
namespace plant {

// The first three state slots of every model, which about fifty readers address
// by position. check_state_layout() is what holds them true.
inline constexpr int HEIGHT_INDEX = 0;
inline constexpr int MORTALITY_INDEX = 1;
inline constexpr int FECUNDITY_INDEX = 2;

// The cumulative hazard at which a cohort is held: dead with probability one,
// its rate parked, and its state still a number.
//
// ⚠️ NEVER PUT A NON-FINITE NUMBER IN AN ODE STATE. `-log(pr_estab)` at an
// establishment probability of exactly zero gives +Inf, and Patch's own
// check_finite_ode_state does not look at this slot -- it tests the cohort
// DENSITY, which is `exp(log_density)` and so a finite zero, and the
// environment's states. Every reader of the hazard is guarded, so a run that
// meets one is correct and silent, and the recorded trajectory it hands the
// reverse sweep carries an entry no derivative can attach to: every census
// metric comes back not-a-number with no refusal declared, which is the one
// failure the gradient must never produce. Measured on a k_I = 20 stand: 9 of
// 88 nodes, from t = 1.5, in 9545 of 11722 recorded steps.
//
// 750 rather than a rounder number because `exp(-x)` is exactly zero past
// 745.14, so survival_individual() and mortality_probability() read the numbers
// +Inf gave them.
//
// ⚠️ EVERY `mortality_dt` MUST TEST THIS BESIDE `is_finite`, and all three do.
// That test is what stops a held cohort's state moving again; keyed on
// finiteness alone a finite hazard passes it, the rate returns, and the run is
// no longer the one +Inf produced.
inline constexpr double establishment_failure_hazard = 750.0;

// The scalar S carries the state, its rates, the auxiliary quantities derived
// from them, and the consumption rates. The consumption rates are five of the
// eleven rate outputs a reverse pass seeds an adjoint on, so storing them at
// their values would read the whole water channel's gradient as zero. The
// environment's water balance is double, and converts where it accumulates.
template <typename S = double>
class Internals {
public:
  using value_type = S;

  Internals(size_t s_size=0, size_t a_size=0, size_t r_size=0)
      :
      states(s_size, S(0.0)),
      rates(s_size, S(NA_REAL)) ,
      auxs(a_size, S(0.0)),
      consumption_rates(r_size, S(NA_REAL))
    {}
  // Perhaps make these private so the () overloads below have some use
  std::vector<S> states;
  std::vector<S> rates;
  std::vector<S> auxs;
  std::vector<S> consumption_rates;  // not quite as pithy

  // Read off the vectors rather than stored beside them. Held as fields they were
  // three more numbers a resize had to remember, and a forgotten one describes a
  // length no vector has.
  size_t state_size() const { return states.size(); }
  size_t aux_size() const { return auxs.size(); }
  size_t resource_size() const { return consumption_rates.size(); }

  // By reference, and at an active scalar that is not a style preference: copying
  // an active scalar registers a tape slot and records an operation, so a read
  // returned by value costs a slot per read per stage for a value nothing writes.
  const S& state(int i) const { return states[i]; }
  const S& rate(int i) const { return rates[i]; }
  const S& aux(int i) const { return auxs[i]; }
  const S& consumption_rate(int i) const { return consumption_rates[i]; }

  void set_state(int i, const S& v) { states[i] = v; }
  void set_rate(int i, const S& v) { rates[i] = v; }
  void set_aux(int i, const S& v) { auxs[i] = v; }
  void set_consumption_rate(int i, const S& v) { consumption_rates[i] = v; }

  void resize(size_t new_size, size_t new_aux_size) {
    states.resize(new_size, S(0.0));
    rates.resize(new_size, S(NA_REAL));
    auxs.resize(new_aux_size, S(0.0));
  }

  void resize_consumption_rates(size_t new_resource_size) {
    consumption_rates.resize(new_resource_size, S(NA_REAL));
  }

  // Every value here that carries the scalar. All four vectors: the rates and
  // the auxiliary quantities as much as the state, because a recording writes
  // them and writing does not refresh what they carry.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, states, rates, auxs, consumption_rates);
  }
};

} // namespace plant

#endif
