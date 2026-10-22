// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <plant/canopy_shape.h> // ShadingModel, shading_model_from_string
#include <odelia/ode_util.hpp>
#include <cmath>                // std::log, std::exp, std::floor (PPA stepping)
#include <type_traits>

using namespace Rcpp;

namespace plant {

// Templated on the scalar S the light field carries; double is production. FF16
// integrates no state of its own, so the field is the whole of what a cohort
// reads from the environment and the whole of what a rebind carries across.
template <typename S = double>
class FF16_Environment : public Environment {
public:
  using value_type = S;

  FF16_Environment() {
    time = 0.0;
  };

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline<S> light_availability;

  // PPA: when true, the light a plant experiences is the stepped (layered)
  // profile rather than the smooth one stored in light_availability. The
  // underlying spline is still fitted to the smooth optical depth (which
  // refines cleanly); the discretisation is applied at read time below.
  bool light_profile_stepped = false;
  // Thickness of one canopy layer in optical-depth units (Control::ppa_layer_optical_depth).
  double layer_optical_depth = 0.5;
  // Smoothing fraction of each layer boundary (Control::ppa_layer_smoothing).
  double layer_smoothing = 0.3;

  // Called once from the Patch constructor. Selects the stepped profile for PPA;
  // deep-crown and crown-centre keep the smooth profile.
  void set_shading_model(const std::string& model,
                         double layer_optical_depth_,
                         double layer_smoothing_) override {
    // Only PPA builds a stepped profile; every other model (including the ""
    // default and mean-light) keeps the smooth profile. Compare the string
    // directly so the "" default does not hit the throwing parser.
    light_profile_stepped = (model == "ppa");
    layer_optical_depth = layer_optical_depth_;
    layer_smoothing = layer_smoothing_;
  }

  // Ability to prescribe a fixed value
  void set_fixed_environment(S value, S height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(S value) {
    S height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  // The smooth profile is returned as read: step_light hands back a copy, and
  // at an active scalar a copy is a tape statement per quadrature point.
  S get_environment_at_height(const S& height) const {
    if (!light_profile_stepped) {
      return light_availability.get_value_at_height(height);
    }
    return step_light(light_availability.get_value_at_height(height));
  }

  // Highest height covered by the light spline; hoist out of hot per-point
  // loops and feed back into the capped get_environment_at_height() overload.
  double max_environment_height() const {
    return light_availability.max_height();
  }

  S get_environment_at_height(const S& height, double cap) const {
    if (!light_profile_stepped) {
      return light_availability.get_value_at_height(height, cap);
    }
    return step_light(light_availability.get_value_at_height(height, cap));
  }

  // Discretise a smooth light value into PPA canopy layers. For the smooth
  // models this is a single predicted branch returning the input unchanged, so
  // it adds no measurable cost to deep-crown/crown-centre. For PPA it maps the
  // optical depth tau = -log(E) onto a smoothed integer number of layers of
  // thickness layer_optical_depth and back-transforms:
  //   E_step = exp(-d * smooth_floor(tau / d)).
  //
  // The stepped profile is double only. Its floor is flat almost everywhere, so
  // a derivative through it would be zero on the flats and a point mass at the
  // joins; a rebind to an active scalar refuses it (assign_from), and this is
  // the backstop for a stepped field that reached an active read anyway.
  S step_light(const S& E) const {
    if (!light_profile_stepped || E >= 1.0) {
      return E;
    }
    // Guard the log: the smooth spline can undershoot to <= 0, which would make
    // tau non-finite. Such a point is fully shaded, so return 0.
    if (E <= 0.0) {
      return S(0.0);
    }
    if constexpr (std::is_same_v<S, double>) {
      const double tau = -std::log(E);
      return std::exp(-layer_optical_depth *
                      smooth_floor(tau / layer_optical_depth));
    } else {
      util::stop("The PPA stepped light profile is a staircase, so it has no "
                 "derivative to carry; take the gradient on a smooth shading "
                 "model");
      return S(0.0); // not reached
    }
  }

  // Monotone, C1-continuous smooth staircase. Each layer is flat over its lower
  // (1 - layer_smoothing) and ramps to the next integer via a cubic smoothstep
  // over its top layer_smoothing fraction. C1 at the joins because smoothstep
  // has zero slope at both ends, so the resulting light profile is smooth enough
  // for the adaptive ODE solver. With layer_smoothing -> 0 this recovers the
  // hard floor (and its instability).
  double smooth_floor(double u) const {
    const double n = std::floor(u);
    const double w = layer_smoothing;
    if (w <= 0.0) {
      return n; // hard step
    }
    const double f = u - n; // fractional position within the layer, [0, 1)
    if (f <= 1.0 - w) {
      return n; // flat lower part of the layer
    }
    const double t = (f - (1.0 - w)) / w; // [0, 1] across the transition
    return n + t * t * (3.0 - 2.0 * t);    // cubic smoothstep
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }
  virtual std::vector<double> get_interpolators_state() const {
    return light_availability.interpolators_state();
  }

  // No integrated state, so no rates of it. A template rather than the base's
  // hidden one, because the patch sums the uptake at the scalar it runs on.
  template <typename V>
  void compute_rates(std::vector<V> const& /* resource_depletion */) {}

  virtual Rcpp::List r_get_state() const {
    return Rcpp::List::create(
              _["light_availability"] = light_availability.r_get_state()
            );
  }

  // The light every knot is left with, from the competition profile above it.
  template <typename Function>
  void compute_environment(Function f_compute_competition_and_slope,
                           S height_max) {
    build_extinction_field(light_availability, f_compute_competition_and_slope,
                           height_max);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }

  // A cohort reads the light field: knot values, then knot slopes, then the
  // canopy top the field is held against. The count and the two fills below
  // are declared together, so a reader that gets the count from this class gets
  // the fill from it too.
  size_t n_cohort_reads() const {
    return 2 * light_availability.knot_count() + 1;
  }

  template <typename It> It cohort_reads(It it) const {
    const std::vector<S>& y = light_availability.knot_values();
    const std::vector<S>& m = light_availability.knot_slopes();
    util::check_length(y.size(), light_availability.knot_count());
    util::check_length(m.size(), light_availability.knot_count());
    for (size_t k = 0; k < y.size(); ++k) {
      util::write_iterator_scalar(it, y[k]);
    }
    for (size_t k = 0; k < m.size(); ++k) {
      util::write_iterator_scalar(it, m[k]);
    }
    util::write_iterator_scalar(it, light_availability.height_max());
    return it;
  }

  template <typename It> It set_cohort_reads(It it) {
    const size_t n_knot = light_availability.knot_count();
    std::vector<S> y(n_knot), m(n_knot);
    for (size_t k = 0; k < n_knot; ++k) { y[k] = *it++; }
    for (size_t k = 0; k < n_knot; ++k) { m[k] = *it++; }
    light_availability.set_knot_data(y, m);
    light_availability.set_height_max(*it++);
    return it;
  }

  // Every active value this environment holds: the light field's knot values
  // and slopes, and the canopy top they are held against.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, light_availability);
  }

  // Another environment's values, written into this one. The light spline is
  // not among them, for the reason a TF24 rebind leaves it: it is derived from
  // the state, and whoever sets a state rebuilds it.
  template <class S1>
  void assign_from(const FF16_Environment<S1>& src) {
    static_cast<Environment&>(*this) = static_cast<const Environment&>(src);
    light_profile_stepped = src.light_profile_stepped;
    layer_optical_depth = src.layer_optical_depth;
    layer_smoothing = src.layer_smoothing;
    if constexpr (!std::is_same_v<S, double>) {
      if (light_profile_stepped) {
        util::stop("The PPA stepped light profile is a staircase, so it has no "
                   "derivative to carry; take the gradient on a smooth shading "
                   "model");
      }
    }
  }

  template <class U>
  FF16_Environment<U> rebind() const {
    FF16_Environment<U> out;
    out.assign_from(*this);
    return out;
  }
};


}

#endif
