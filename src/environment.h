// -*-c++-*-
#ifndef TREE_ENVIRONMENT_H_
#define TREE_ENVIRONMENT_H_

#include "parameters.h"
#include "spline.h"
#include "util.h"
#include "functor.h"

namespace model {

class Environment {
public:
  Environment(Parameters p);
  double canopy_openness(double height) const;
  void compute_light_environment(util::DFunctor *canopy_openness,
				 double height_max);
  double patch_survival(double time_at_birth) const;
  void clear();

  // NOTE: Interface here will change
  double seed_rain_rate() const;
  void set_seed_rain_index(size_t x);

  double get_time() const;
  void set_time(double x);

  // * R interface
  std::vector<double> r_get_seed_rain() const;
  void r_set_seed_rain(std::vector<double> x);

  spline::Spline r_get_light_environment() const;
  void r_set_light_environment(const spline::Spline env);

  void r_set_seed_rain_index(size_t x);

private:
  spline::Spline light_environment;
  Disturbance disturbance_regime;
  std::vector<double> seed_rain;
  size_t seed_rain_index;
  Control control;
  double time;
};

}

RCPP_EXPOSED_CLASS(model::Environment)

#endif
