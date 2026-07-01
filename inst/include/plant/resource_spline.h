// -*-c++-*-

#ifndef PLANT_PLANT_RESOURCE_SPLINE_H_
#define PLANT_PLANT_RESOURCE_SPLINE_H_

#include <odelia/interpolator.hpp>
#include <plant/adaptive_interpolator.h>
#include <odelia/ode_interface.hpp>
#include <plant/util.h>
#include <algorithm> // std::max, for the resource-availability floor (#253)

using namespace Rcpp;

namespace plant {

class ResourceSpline {
public:

  // Constructors
  ResourceSpline() {
    setup(1e-6, 17, 16, false);
  }

  ResourceSpline(double tol, size_t nbase, size_t max_depth, bool rescale_usually) {
    setup(tol, nbase, max_depth, rescale_usually);
  }

  void setup(double tol, size_t nbase, size_t max_depth, bool rescale_usually) {

    // Initialise adaptive interpolator. This object can create an interpolator spline
    spline_construction =
        interpolator::AdaptiveInterpolator(tol, tol, nbase, max_depth);
    // Create an actual spline, For initalisation
    // Provide a dummy function and construct
    // This will be over-written later with actual function
    spline = spline_construction.construct(
        [&](double height) { return get_value_at_height(height); }, 0,
        1); // these are update with init(x, y) when patch is created
    
    spline_rescale_usually = rescale_usually;
  };

  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {
    if (rescale & spline_rescale_usually) {
      rescale_spline(f_compute_competition, height_max);
    } else {
      construct_spline(f_compute_competition, height_max);
    }
  };

  void set_fixed_value(double value, double height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<double> y = {value, value, value};
    clear();
    spline.init(x, y);
  }

  void clear() {
    spline.clear();
  }

  // Highest height covered by the spline; above this get_value_at_height()
  // returns the hard-coded open value (1.0). Hoist this out of hot per-point
  // loops with get_value_at_height(height, cap).
  double max_height() const { return spline.max(); }

  double get_value_at_height(double height) const {
    return get_value_at_height(height, spline.max());
  }

  // Variant taking a pre-fetched cap (= max_height()) so callers integrating
  // over many points pay the spline.max() lookup once rather than per point.
  double get_value_at_height(double height, double cap) const {
    // TODO(#385): change maximum - here hard-coded to 1.0
    // `cap` already guards the upper bound and the crown integral keeps
    // height >= 0 = spline.min(), so use the unchecked operator() rather than
    // eval() to avoid re-running check_active()/bound checks per quadrature
    // point. Same underlying tk_spline(height) call.
    //
    // Floor the result at 0 (#253): the cubic spline can undershoot below
    // zero between knots (notably the K93 light spline at high k_I), which is
    // non-physical for a resource availability. The clamp is a no-op for the
    // usual positive case (so values stay bit-identical there) and only bites
    // on the spurious negative undershoot. This is the single chokepoint for
    // FF16/K93/TF24, and belongs here in plant rather than in the
    // general-purpose interpolator (which is migrating to odelia).
    return height <= cap ? std::max(0.0, spline(height)) : 1.0;
  }

  // Analytic d(value)/d(height), consistent with get_value_at_height: 0 above
  // the cap (constant open value) and where the #253 floor clamps a negative
  // undershoot, else the spline's analytic derivative. Enables exact AD
  // gradients through the resource environment (#472 scope B / #537).
  double get_value_deriv_at_height(double height) const {
    return get_value_deriv_at_height(height, spline.max());
  }
  double get_value_deriv_at_height(double height, double cap) const {
    if (height > cap) return 0.0;
    return spline(height) > 0.0 ? spline.deriv(height) : 0.0;
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

  // This object will store an interpolator spline of 
  // resource availability as a function of size
  odelia::interpolator::Interpolator spline;

  // This object can create an interpolator spline via adaptive refinement
  interpolator::AdaptiveInterpolator spline_construction;

  // flag, do we try to rescale the spline when possible? this is quicker
  bool spline_rescale_usually;

  Rcpp::NumericMatrix r_get_state() const
  {

    // format spline as Matrix
    std::vector<std::vector<double>> xy;
    xy.push_back(spline.get_x());
    xy.push_back(spline.get_y());

    const size_t n = xy.size();
    Rcpp::NumericMatrix ret(static_cast<int>(xy.begin()->size()),
                            static_cast<int>(n));
    Rcpp::NumericMatrix::iterator it = ret.begin();
    for (size_t i = 0; i < n; ++i)
    {
      it = std::copy(xy[i].begin(), xy[i].end(), it);
    }

    // Add colnames
    ret.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::CharacterVector::create("height", "light_availability"));
    return ret;
  }

private:

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

  };


} // plant namespace

#endif
