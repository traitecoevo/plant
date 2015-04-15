// -*-c++-*-
#ifndef TREE_ENVIRONMENT_H_
#define TREE_ENVIRONMENT_H_

#include "parameters.h"
#include "interpolator.h"
#include "util.h"
#include "functor.h"
#include "adaptive_interpolator.h"

namespace model {

class Environment {
public:
  Environment(Parameters p);
  double canopy_openness(double height) const;
  void compute_light_environment(util::DFunctor *f_canopy_openness,
				 double height_max);
  void rescale_light_environment(util::DFunctor *f_canopy_openness,
				 double height_max);
  double patch_survival() const;
  double patch_survival_conditional(double time_at_birth) const;
  void clear();

  // NOTE: Interface here will change
  double seed_rain_rate() const;
  void set_seed_rain_index(size_t x);

  const Disturbance& get_disturbance_regime() const;

  double get_time() const;
  void set_time(double x);

  interpolator::Interpolator get_light_environment() const;
  void set_light_environment(const interpolator::Interpolator env);

  // * R interface
  std::vector<double> r_get_seed_rain() const;
  void r_set_seed_rain(std::vector<double> x);

  void r_set_seed_rain_index(size_t x);

  Rcpp::List r_get_state() const;
  void r_set_state(Rcpp::List x);

private:
  interpolator::Interpolator light_environment;
  Disturbance disturbance_regime;
  std::vector<double> seed_rain;
  size_t seed_rain_index;
  Control control;
  double time;
  interpolator::AdaptiveInterpolator light_environment_generator;
};

}

RCPP_EXPOSED_CLASS_NODECL(model::Environment)

#endif
