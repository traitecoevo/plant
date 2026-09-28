// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_

#include <plant/control.h>
#include <odelia/ode_interface.hpp>
#include <plant/internals.h>
#include <plant/util.h>
#include <unordered_map>
#include <Rcpp.h>
#include <plant/extrinsic_drivers.h>

using namespace Rcpp;

namespace plant {

class Environment {
public:
  // Configure the crown shading model for the light profile. Default: no-op;
  // only FF16_Environment builds an alternative (stepped) profile. Called once
  // from the Patch constructor with the run's Control settings.
  virtual void set_shading_model(const std::string& /*model*/,
                                 double /*layer_optical_depth*/,
                                 double /*layer_smoothing*/) {}

  // ODE interface. An environment holding no integrated state answers zero and
  // leaves every iterator where it found it. One that integrates state declares
  // these itself, at the scalar it holds that state in, so the state and the
  // rates a cohort reads back travel at the same scalar the cohort does.
  size_t ode_size() const { return 0; }

  // How many entries of resource_depletion compute_rates reads, and so how
  // long each individual's consumption vector must be.
  virtual size_t n_resources() const { return 0; }

  // One aux slot per resource, holding the uptake the individuals supplied to
  // compute_rates. The soil rates subtract it and keep no record of it, so
  // without this slot the consumption is unrecoverable from the state.
  size_t aux_size() const { return n_resources(); }

  // Add `amount` of resource `i` at one instant (#628). What the resource is,
  // and what the amount is measured in, is the environment's business: TF24's
  // resources are soil layers and the amount is metres of water.
  //
  // Applied between solver legs, so nothing about it is error-controlled -- an
  // implementation must bound the jump itself, because no error estimate and
  // no step rejection stand behind it. Returns what it managed to apply:
  // {accepted, shed}, in the same units.
  //
  // An environment with no resources refuses rather than silently swallowing
  // the amount, so a pulse aimed at a model that cannot take one is visible.
  virtual std::vector<double> add_resource_pulse(size_t /*i*/,
                                                 double /*amount*/) {
    util::stop("This environment has no resource pools, so there is nothing "
               "for a resource pulse to add to");
    return std::vector<double>(); // not reached
  }

  template <typename It> It set_ode_state(It it) { return it; }

  template <typename It> It ode_state(It it) const { return it; }

  template <typename It> It ode_rates(It it) const { return it; }

  template <typename It> It ode_aux(It it) const { return it; }

  template <typename It> It set_ode_aux(It it) { return it; }

  // No integrated state, so no rates of it. An environment that has them
  // declares its own, taking the uptake at the scalar the patch summed it in.
  template <typename V> void compute_rates(const std::vector<V>&) {}

  virtual Rcpp::List r_get_state() const
  {
    return Rcpp::List::create(_["time"] = time);
  }

  // Reset the environment
  void clear() {
    time = 0.0;
    clear_environment();
    clear_state();
  }

  // Discard the competition profile. A patch with no individuals casts no
  // shade, so this runs whenever the patch empties, not only between runs.
  virtual void clear_environment() {}

  // Restore integrated state to the values a run starts from. Only clear()
  // calls this: an empty patch has nothing to shade with but the state the
  // solver has integrated is still live.
  //
  // These are two virtuals rather than one because they have different
  // lifetimes, not because any environment needs both. StochasticPatch clears
  // the profile on an empty patch; the state has to survive that. Note that
  // Patch does *not* clear -- it skips recomputation and leaves the profile
  // standing -- so the deterministic path never exercises the distinction, and
  // making StochasticPatch behave like Patch would have avoided the split.
  // Keeping it means an environment that acquires ODE state later cannot get
  // this wrong by inheriting a clear_environment() that resets it.
  virtual void clear_state() {}

  virtual void r_init_interpolators(const std::vector<double>& state) {}
  // What r_init_interpolators() reads back. An environment with no interpolant
  // has none, which is what makes this a no-op pair rather than a special case.
  virtual std::vector<double> get_interpolators_state() const { return {}; }

  double get_environment_at_height(double height) const { return 0.0; };

  virtual ~Environment() = default;

  double time;

  size_t species_arriving_index;

  ExtrinsicDrivers extrinsic_drivers;

  // The
  std::vector<std::string> extrinsic_drivers_get_names() const
  {
    return  extrinsic_drivers.get_names();
  }

  void extrinsic_drivers_set_constant(std::string driver_name, double value)
  {
    extrinsic_drivers.set_constant(driver_name, value);
  }

  // ⚠️ MONOTONE, NOT odelia's natural default. A driver is data given as values
  // alone, and a natural spline beside a wet day reads NEGATIVE rainfall between
  // the points (test-tf24-water-budget.R); monotone slopes keep every read inside
  // the two values bracketing it. The birth-rate drivers make the same choice.
  void extrinsic_drivers_set_variable(std::string driver_name, std::vector<double> const &x, std::vector<double> const &y)
  {
    extrinsic_drivers.set_variable(driver_name, x, y, odelia::drivers::Slopes::monotone);
  }

  double extrinsic_drivers_evaluate(std::string driver_name, double x) const
  {
    return extrinsic_drivers.evaluate(driver_name, x);
  }

  std::vector<double> extrinsic_drivers_evaluate_range(std::string driver_name, std::vector<double> const &x) const
  {
    return extrinsic_drivers.evaluate_range(driver_name, x);
  }
  
};
}
#endif
