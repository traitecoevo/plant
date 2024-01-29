// -*-c++-*-
// todo: What does this class do?

#ifndef PLANT_PLANT_SHADING_H_
#define PLANT_PLANT_SHADING_H_

#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/ode_interface.h>
#include <plant/util.h>

using namespace Rcpp;

namespace plant {

class Shading {
public:

  // Constructors
  Shading() {
    
    // Initialise adaptive interpolator. This object can create an interpolator spline
    spline_construction = interpolator::AdaptiveInterpolator(1e-6, 1e-6, 17, 16);
    // Create an actual spline, For initalisation
    // Provide a dummy function and construct
    // This will be over-written later with actual function
    spline = spline_construction.construct(
      [&](double height) {
          return get_canopy_at_height(height);
      }, 0, 1); // these are update with init(x, y) when patch is created

    shading_spline_rescale_usually = false;
  };

  Shading(double tol, size_t nbase, size_t max_depth, bool rescale_usually) {

    // Initialise adaptive interpolator. This object can create an interpolator spline
    spline_construction =
        interpolator::AdaptiveInterpolator(tol, tol, nbase, max_depth);
    // Create an actual spline, For initalisation
    // Provide a dummy function and construct
    // This will be over-written later with actual function
    spline = spline_construction.construct(
        [&](double height) { return get_canopy_at_height(height); }, 0,
        1); // these are update with init(x, y) when patch is created
    
    shading_spline_rescale_usually = rescale_usually;
  };

  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {
    if (rescale & shading_spline_rescale_usually) {
      rescale_spline(f_compute_competition, height_max);
    } else {
      construct_spline(f_compute_competition, height_max);
    }
  };

  template <typename Function>
  void construct_spline(Function f_compute_competition, double height_max)
  {
    const double lower_bound = 0.0;
    double upper_bound = height_max;

    spline =
      spline_construction.construct(f_compute_competition, lower_bound, upper_bound);
  }

  template <typename Function>
  void rescale_spline(Function f_compute_competition, double height_max) {
    std::vector<double> h = spline.get_x();
    const double min = spline.min(), // 0.0?
      height_max_old = spline.max();

    util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
    h.back() = height_max; // Avoid round-off error.

    spline.clear();
    for (auto hi : h) {
      spline.add_point(hi, f_compute_competition(hi));
    }
    spline.initialise();
  }

  void set_fixed_values(double value, double height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<double> y = {value, value, value};
    clear();
    spline.init(x, y);
  }

  void clear() {
    spline.clear();
  }

  double get_canopy_at_height(double height) const {
    const bool within = height <= spline.max();
    // TODO: change maximum - here hard-coded to 1.0
    return within ? spline.eval(height) : 1.0;
  }

  virtual void r_init_interpolators(const std::vector<double>& state) {
    // See issue #144; this is important as we have to at least refine
    // the light environment, but doing this is better because it means
    // that if rescale_usually is on we do get the same light
    // environment as before.
    if (state.size() % 2 != 0) {
      util::stop("Expected even number of elements in light environment");
    }
    const size_t state_n = state.size() / 2;
    auto it = state.begin();
    std::vector<double> state_x, state_y;
    std::copy_n(it,         state_n, std::back_inserter(state_x));
    std::copy_n(it + state_n, state_n, std::back_inserter(state_y));
    spline.init(state_x, state_y);
  }

    // This object will store an interpolator spline of shading
    interpolator::Interpolator spline;

    // This object can create an interpolator spline via adaptive refinement
    interpolator::AdaptiveInterpolator spline_construction;

    // flag, do we try to rescale the spline when possible? this is quicker
    bool shading_spline_rescale_usually;
  };

inline Rcpp::NumericMatrix get_state(const Shading shading) {
  using namespace Rcpp;
  NumericMatrix xy = shading.spline.r_get_xy();
  Rcpp::CharacterVector colnames =
    Rcpp::CharacterVector::create("height", "canopy_openness");
  xy.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);
  return xy;
}


}

#endif
