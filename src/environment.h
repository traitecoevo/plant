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
  double patch_survival(double age_at_birth) const;
  void clear();

  spline::Spline get_light_environment() const;
  void set_light_environment(const spline::Spline env);

  double get_age() const;
  void set_age(double x);

private:
  spline::Spline light_environment;
  Disturbance disturbance_regime;
  Control control;
  double age;
};

}

RCPP_EXPOSED_CLASS(model::Environment)

#endif
