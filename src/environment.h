// -*-c++-*-
#ifndef TREE_ENVIRONMENT_H_
#define TREE_ENVIRONMENT_H_

#include "spline.h"
#include "disturbance.h"
#include "util.h"
#include "functor.h"
#include "control.h"

namespace model {

class Environment {
public:
  // TODO: Allow (require?) initialisation from Parameters
  // TODO: clear method (invoke from Patch::clear)
  Environment(Disturbance disturbance_regime,
	      Control control);
  double canopy_openness(double height) const;
  void compute_light_environment(util::DFunctor *canopy_openness,
				 double height_max);
  double patch_survival(double age_at_birth) const;

  spline::Spline get_light_environment() const;
  void set_light_environment(const spline::Spline env);

  Disturbance get_disturbance_regime() const;

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
