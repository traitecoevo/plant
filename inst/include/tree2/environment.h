// -*-c++-*-
#ifndef TREE_ENVIRONMENT_H_
#define TREE_ENVIRONMENT_H_

#include <tree2/control.h>
#include <tree2/disturbance.h>
#include <tree2/interpolator.h>
#include <tree2/adaptive_interpolator.h>
#include <tree2/util.h>

namespace tree2 {

class Environment {
public:
  template <typename Parameters>
  Environment(Parameters p);
  double canopy_openness(double height) const;
  template <typename Function>
  void compute_light_environment(Function f_canopy_openness, double height_max);
  template <typename Function>
  void rescale_light_environment(Function f_canopy_openness, double height_max);
  double patch_survival() const;
  double patch_survival_conditional(double time_at_birth) const;
  void clear();

  // NOTE: Interface here will change
  double seed_rain_dt() const;
  void set_seed_rain_index(size_t x);

  // * R interface
  void r_set_seed_rain_index(util::index x);

  double time;
  Disturbance disturbance_regime;
  interpolator::Interpolator light_environment;

private:
  std::vector<double> seed_rain;
  size_t seed_rain_index;
  interpolator::AdaptiveInterpolator light_environment_generator;
};

template <typename Parameters>
Environment::Environment(Parameters p)
  : time(0.0),
    disturbance_regime(p.disturbance_mean_interval),
    seed_rain(p.seed_rain),
    seed_rain_index(0),
    light_environment_generator(make_interpolator(p.control)) {
}

template <typename Function>
void Environment::compute_light_environment(Function f_canopy_openness,
                                            double height_max) {
  light_environment =
    light_environment_generator.construct(f_canopy_openness, 0, height_max);
}

template <typename Function>
void Environment::rescale_light_environment(Function f_canopy_openness,
                                            double height_max) {
  std::vector<double> h = light_environment.get_x();
  const double min = light_environment.min(), // 0.0?
    height_max_old = light_environment.max();

  util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
  h.back() = height_max; // Avoid round-off error.

  light_environment.clear();
  for (auto hi : h) {
    light_environment.add_point(hi, f_canopy_openness(hi));
  }
  light_environment.initialise();
}

inline interpolator::AdaptiveInterpolator
make_interpolator(const Control& control) {
  using namespace interpolator;
  return AdaptiveInterpolator(control.environment_light_tol,
                              control.environment_light_tol,
                              control.environment_light_nbase,
                              control.environment_light_max_depth);
}

}

#endif
