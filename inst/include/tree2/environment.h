// -*-c++-*-
#ifndef TREE_ENVIRONMENT_H_
#define TREE_ENVIRONMENT_H_

#include <tree2/parameters.h>
#include <tree2/interpolator.h>
#include <tree2/adaptive_interpolator.h>
#include <tree2/util.h>

namespace tree2 {

class Environment {
public:
  Environment(Parameters p);
  double canopy_openness(double height) const;
  // void compute_light_environment(util::DFunctor *f_canopy_openness,
  //                                double height_max);
  // void rescale_light_environment(util::DFunctor *f_canopy_openness,
  //                                double height_max);
  double patch_survival() const;
  double patch_survival_conditional(double time_at_birth) const;
  void clear();

  // NOTE: Interface here will change
  double seed_rain_rate() const;
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

}

#endif
