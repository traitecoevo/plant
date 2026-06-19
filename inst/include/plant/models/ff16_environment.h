// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <plant/interpolator.h>
#include <plant/canopy_shape.h> // ShadingModel, shading_model_from_string
#include <cmath>                // std::log, std::exp, std::floor (PPA stepping)

using namespace Rcpp;

namespace plant {

class FF16_Environment : public Environment {
public:
  // constructor for R interface - default settings can be modified
  // except for light_availability_spline_rescale_usually
  // which are only updated on construction
  FF16_Environment() {
    time = 0.0;

    // Shading defaults have lower tolerance which are overwritten for speed
    light_availability = ResourceSpline(
        1e-4, // light_availability_spline_tol,
        17,   // light_availability_spline_nbase,
        16,   // light_availability_spline_max_depth,
        true  // light_availability_spline_rescale_usually)
    );

  };

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline light_availability;

  // PPA: when true, the light a plant experiences is the stepped (layered)
  // profile rather than the smooth one stored in light_availability. The
  // underlying spline is still fitted to the smooth optical depth (which
  // refines cleanly); the discretisation is applied at read time below.
  bool light_profile_stepped = false;
  // Thickness of one canopy layer in optical-depth units (Control::ppa_layer_optical_depth).
  double layer_optical_depth = 0.5;

  // Called once from the Patch constructor. Selects the stepped profile for PPA;
  // deep-crown and flat-top keep the smooth profile.
  void set_shading_model(const std::string& model,
                         double layer_optical_depth_) override {
    light_profile_stepped = (shading_model_from_string(model) == ShadingModel::PPA);
    layer_optical_depth = layer_optical_depth_;
  }

  // Ability to prescribe a fixed value
  // TODO(#476): add setting to set other variables like water
  void set_fixed_environment(double value, double height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return step_light(light_availability.get_value_at_height(height));
  }

  // Highest height covered by the light spline; hoist out of hot per-point
  // loops and feed back into the capped get_environment_at_height() overload.
  double max_environment_height() const {
    return light_availability.max_height();
  }

  double get_environment_at_height(double height, double cap) const {
    return step_light(light_availability.get_value_at_height(height, cap));
  }

  // Discretise a smooth light value into PPA canopy layers. For the smooth
  // models this is a single predicted branch returning the input unchanged, so
  // it adds no measurable cost to deep-crown/flat-top. For PPA it floors the
  // optical depth tau = -log(E) to an integer number of layers of thickness
  // layer_optical_depth and back-transforms: E_step = exp(-d * floor(tau / d)).
  double step_light(double E) const {
    if (!light_profile_stepped || E >= 1.0) {
      return E;
    }
    // Guard the log: the smooth spline can undershoot to <= 0, which would make
    // tau non-finite. Such a point is fully shaded, so return 0.
    if (E <= 0.0) {
      return 0.0;
    }
    const double tau = -std::log(E);
    return std::exp(-layer_optical_depth * std::floor(tau / layer_optical_depth));
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }

  virtual void compute_rates(std::vector<double> const& resource_depletion) {

  }

  virtual Rcpp::List r_get_state() const {
    return Rcpp::List::create(
              _["light_availability"] = light_availability.r_get_state()
            );
  }

  // Pre-compute resources available in the environment, as a function of height
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of light_availability spline
    // Note: extinction coefficient was already applied in strategy, so
    // f_compute_competition gives sum of projected leaf area (k L) across species. Just need to apply Beer's law, E = exp(- (k L))
    auto f_light_availability = [&](double height) -> double
    { return exp(-f_compute_competition(height)); };

    // Calculates the light_availability spline, by fitting to the function
    // `f_compute_competition` as a function of height
    light_availability.compute_environment(f_light_availability, height_max, rescale);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }
};


}

#endif
